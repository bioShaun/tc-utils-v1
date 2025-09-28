import io
import zipfile
from pathlib import Path
from typing import Optional

import polars as pl
import typer
from loguru import logger


def load_one_month_data(zip_filename: Path) -> pl.DataFrame:
    """Loads data for one month from a zip file containing Excel files, using Polars."""
    all_data_frames = []

    try:
        # Change extension to .ipc for Polars IPC format
        zip_filename_ipc = zip_filename.with_suffix(".ipc")
        # Check if an IPC file exists
        if zip_filename_ipc.exists():
            logger.info(f"{zip_filename_ipc} exists, loading from IPC file")
            # Read the IPC file and select necessary columns
            return pl.read_ipc(zip_filename_ipc).select(
                ["扣费时间", "退出状态", "CPU核*时", "金额(元)", "开始时间", "结束时间"]
            )

        # Check if the zip file exists
        if not zip_filename.exists():
            raise FileNotFoundError(f"错误：找不到 Zip 文件 '{zip_filename}'。")

        with zipfile.ZipFile(zip_filename, "r") as zf:
            # Find all .xlsx files in the zip archive
            excel_files = [f for f in zf.namelist() if f.lower().endswith(".xlsx")]

            if not excel_files:
                logger.warning(f"错误：在 {zip_filename} 中没有找到 .xlsx 文件。")
            else:
                logger.info(f"找到 Excel 文件: {excel_files}")
                # Read each Excel file
                for excel_file_name in excel_files:
                    try:
                        logger.info(f"正在读取 '{excel_file_name}'...")
                        # Open the Excel file from the zip archive
                        with zf.open(excel_file_name) as f:
                            # Read file content into an in-memory bytes stream
                            excel_content = io.BytesIO(f.read())
                            # Use Polars to read the Excel file (might need 'xlsx2csv' installed: pip install polars[excel])
                            # Read all sheets by setting sheet_id=0 (0 means all sheets)
                            # sheet_id=0 returns a dict of {sheet_name: DataFrame}
                            sheet_data_dict = pl.read_excel(excel_content, sheet_id=0)
                            # Extend the list with all DataFrames from the sheets
                            all_data_frames.extend(sheet_data_dict.values())
                        logger.info(f"完成读取 '{excel_file_name}'。")
                    except Exception as e:
                        logger.error(f"读取文件 '{excel_file_name}' 时出错: {e}")

        # If data was successfully read
        if all_data_frames:
            # Concatenate data from all sheets and files
            final_df = pl.concat(all_data_frames)
            logger.info(f"将工作表存储到 {zip_filename_ipc}")
            # Save the combined DataFrame to an IPC file
            final_df.write_ipc(zip_filename_ipc)
            # Select and return the required columns
            return final_df.select(
                ["扣费时间", "退出状态", "CPU核*时", "金额(元)", "开始时间", "结束时间"]
            )
        logger.warning("\n未能读取任何数据，未生成输出文件。")

    except FileNotFoundError as e:
        logger.error(e)
    except zipfile.BadZipFile:
        logger.error(f"错误：'{zip_filename}' 不是有效的 zip 文件或已损坏。")
    except ImportError:
        logger.error("错误：缺少必要的库 (polars, xlsx2csv)。")
        logger.error("请先安装它们： pip install polars[excel]")
    except Exception as e:
        logger.error(f"发生意外错误: {e}")
    # Return an empty DataFrame if errors occur
    return pl.DataFrame()


def one_month_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Calculates monthly statistics using Polars."""
    if df.is_empty():
        logger.warning("没有数据可供处理。")
        return pl.DataFrame()

    # Convert '扣费时间' column to datetime, handling potential errors
    # Assuming the format might be like 'YYYY-MM-DD HH:MM:SS' or similar
    # Polars might infer the type, but explicit conversion is safer.
    # Use try_parse_dates for robustness if formats vary.
    df = df.with_columns(
        pl.col("扣费时间").str.to_datetime(
            strict=False
        )  # Use strict=False to allow parsing various formats or nulls
    )

    # Filter out rows with invalid dates (null after conversion)
    df = df.filter(pl.col("扣费时间").is_not_null())

    # Extract year-month information and format as string 'YYYY-MM'
    df = df.with_columns(pl.col("扣费时间").dt.strftime("%Y-%m").alias("年月"))

    # Explicitly cast aggregation columns to Float64 before summing
    # Use strict=False to turn conversion errors into nulls instead of raising an error
    # Fill any nulls resulting from the cast with 0 before summing
    df = df.with_columns(
        pl.col("CPU核*时").cast(pl.Float64, strict=False).fill_null(0),
        pl.col("金额(元)").cast(pl.Float64, strict=False).fill_null(0),
    )

    # Add logging to inspect DataFrame before aggregation
    logger.debug(f"DataFrame shape before aggregation: {df.shape}")
    logger.debug(f"DataFrame dtypes before aggregation: {df.dtypes}")
    # Optional: Log head for more detail, might be verbose
    # logger.debug(f"DataFrame head before aggregation:\n{df.head()}")

    # Calculate monthly total CPU core hours and total amount per '退出状态'
    monthly_stats = df.group_by(["年月", "退出状态"]).agg(
        [
            pl.sum("CPU核*时").alias("CPU核*时"),
            pl.sum("金额(元)").alias("金额(元)"),
        ]
    )

    return monthly_stats


def main(stats_dir: Path, summary_filename: Path, prefix: Optional[str] = None):
    """Main function to process all zip files in a directory and generate a summary CSV."""
    if not stats_dir.exists():
        logger.error(f"错误：找不到目录 '{stats_dir}'。")
        return

    stats_df_list = []
    # Iterate over all .zip files in the specified directory
    if prefix is None:
        stats_files = stats_dir.glob("*.zip")
    else:
        stats_files = stats_dir.glob(f"*{prefix}*zip")
    for file_i in stats_files:
        logger.info(f"正在处理文件 {file_i}")
        df = load_one_month_data(file_i)
        if not df.is_empty():
            stats_df = one_month_stats(df)
            stats_df_list.append(stats_df)

    if not stats_df_list:
        logger.warning("没有找到或处理任何有效的月份数据。")
        return

    logger.info("正在合并所有月份的统计数据")
    # Concatenate statistics from all months
    final_stats_df = pl.concat(stats_df_list)

    # Explicitly cast 'CPU核*时' to Float64 before calculation
    final_stats_df = final_stats_df.with_columns(
        pl.col("CPU核*时").cast(
            pl.Float64, strict=False
        )  # Use strict=False to handle potential conversion errors gracefully (produces null)
    )

    # Calculate the '实际金额(元)' based on 'CPU核*时'
    final_stats_df = final_stats_df.with_columns(
        (pl.col("CPU核*时") * 0.04).alias("实际金额(元)")
    )

    # Save the final statistics to a CSV file
    logger.info(f"正在将最终统计数据写入 {summary_filename}")
    # final_stats_df.write_csv(summary_filename, separator="\t")
    out_final_stats_df = final_stats_df.to_pandas()
    out_final_stats_df.to_excel(summary_filename, index=False)
    logger.info("处理完成。")


if __name__ == "__main__":
    typer.run(main)
