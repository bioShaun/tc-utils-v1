from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import pandas as pd
import pandera as pa
from loguru import logger
from pandera.typing import Series

PROBE_COLUMNS = ["chrom", "probe_start", "probe_end", "pos"]


class DesignDBSchema(pa.DataFrameModel):
    chrom: Series[str] = pa.Field(nullable=False)
    pos: Series[int] = pa.Field(nullable=False)
    probe_start: Series[int] = pa.Field(nullable=False)
    probe_end: Series[int] = pa.Field(nullable=False)


@dataclass
class AnnDataProcessor:
    ann_table: Path
    validated_df: pd.DataFrame = field(init=False)

    def load_data(self, schema: pa.DataFrameModel) -> None:
        """加载并验证数据"""
        try:
            df = pd.read_table(
                self.ann_table,
            )
            self.validated_df = schema.validate(df)  # type: ignore
            logger.info("数据加载和验证成功")
        except pa.errors.SchemaError as e:  # type: ignore
            logger.error(f"数据验证失败: {str(e)}")
            raise
        except Exception as e:
            logger.error(f"数据加载失败: {str(e)}")
            raise

    @lru_cache(maxsize=1)
    def get_validated_data(self) -> pd.DataFrame:
        """缓存验证后的数据"""
        return self.validated_df

    def to_probe_bed(self):
        """将ann_table转换为probe_bed"""
        probe_bed = self.ann_table.with_suffix(".panel.bed")
        self.validated_df.sort_values(["chrom", "probe_start"], inplace=True)
        self.validated_df.to_csv(
            probe_bed, sep="\t", index=False, header=False, columns=PROBE_COLUMNS
        )
        logger.info(f"生成探针BED文件: {probe_bed}")
        return probe_bed

    def to_target_bed(self):
        pass

    def to_target_id(self):
        pass
