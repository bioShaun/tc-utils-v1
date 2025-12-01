from pathlib import Path
from typing import Tuple
import pandas as pd
import typer
from loguru import logger

def file2excel(input_path: Path, excel_path: Path) -> Tuple[int, int]:
    """
    Convert a CSV or TSV file to an Excel file.

    Args:
        input_path (Path): Path to the input CSV or TSV file.
        excel_path (Path): Path to the output Excel file.

    Returns:
        Tuple[int, int]: The shape of the dataframe (rows, columns).

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the file format cannot be determined or parsed.
    """
    if not input_path.exists():
        logger.error(f"File not found: {input_path}")
        raise FileNotFoundError(f"File not found: {input_path}")

    logger.info(f"Reading file: {input_path}")
    
    try:
        # Determine separator based on file extension
        if input_path.suffix.lower() == '.tsv':
            df = pd.read_table(input_path)
        elif input_path.suffix.lower() == '.csv':
            df = pd.read_csv(input_path)
        else:
            # Attempt to auto-detect separator
            logger.debug("Extension not recognized, attempting to auto-detect separator.")
            try:
                with open(input_path, 'r', encoding='utf-8') as f:
                    first_line = f.readline()
                
                if '\t' in first_line:
                    df = pd.read_table(input_path)
                elif ',' in first_line:
                    df = pd.read_csv(input_path)
                else:
                    # Fallback to python engine with auto detection
                    df = pd.read_csv(input_path, sep=None, engine='python')
            except Exception as e:
                logger.error(f"Failed to auto-detect format: {e}")
                raise ValueError(f"Could not determine file format for {input_path}") from e

        logger.info(f"Data loaded. Shape: {df.shape}")
        
        # Save to Excel
        df.to_excel(excel_path, index=False)
        logger.info(f"Successfully saved to {excel_path}")
        
        return df.shape

    except Exception as e:
        logger.error(f"An error occurred during conversion: {e}")
        raise

if __name__ == "__main__":
    typer.run(file2excel)
