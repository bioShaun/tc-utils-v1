from pathlib import Path
from typing import Annotated

import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta
from rich.console import Console
from tqdm import tqdm

console = Console()


def validate_file(path: Path, name: str) -> None:
    """Validate that a path exists and is a file."""
    if not path.exists():
        logger.error(f"{name} file not found: {path}")
        console.print(f"[red]Error:[/red] {name} file not found: {path}")
        raise typer.Exit(code=1)

    if not path.is_file():
        logger.error(f"{name} path is not a file: {path}")
        console.print(f"[red]Error:[/red] {name} path is not a file: {path}")
        raise typer.Exit(code=1)


def load_id_set(id_file: Path) -> set[str]:
    """Load IDs from a text file, one ID per line."""
    validate_file(id_file, "ID list")
    id_set = {line.strip() for line in id_file.read_text(encoding="utf-8").splitlines() if line.strip()}

    if not id_set:
        logger.error(f"ID list file is empty: {id_file}")
        console.print(f"[red]Error:[/red] ID list file is empty: {id_file}")
        raise typer.Exit(code=1)

    logger.info(f"Loaded {len(id_set)} IDs from {id_file}")
    return id_set


def main(
    vcf_file: Annotated[Path, typer.Argument(help="Input VCF file path")],
    ref: Annotated[Path, typer.Argument(help="Reference FASTA file path")],
    out_file: Annotated[Path, typer.Argument(help="Output TSV file path")],
    flank_size: Annotated[int, typer.Option("--flank-size", help="Flanking sequence size")] = 200,
    id_file: Annotated[
        Path | None,
        typer.Option("--id", help="Optional ID list file path; one ID per line"),
    ] = None,
) -> None:
    """Generate primer sequences from VCF and reference files."""
    validate_file(vcf_file, "VCF")
    validate_file(ref, "Reference")
    id_set: set[str] | None = load_id_set(id_file) if id_file else None

    logger.info("Reading VCF file...")
    vcf_df = pd.read_table(
        vcf_file,
        comment="#",
        usecols=[0, 1, 2, 3, 4],
        names=["chrom", "pos", "id", "ref", "alt"],
    )
    vcf_df["chrom"] = vcf_df["chrom"].astype("str")
    vcf_df["id"] = vcf_df["id"].astype("str")

    if (vcf_df["id"] == ".").all():
        vcf_df["effective_id"] = vcf_df["chrom"] + "_" + vcf_df["pos"].astype("str")
    else:
        vcf_df["effective_id"] = vcf_df["id"]

    if id_set is not None:
        vcf_df = vcf_df[vcf_df["effective_id"].isin(id_set)]
        logger.info(f"Filtered variants by ID list, {len(vcf_df)} variants remain")
        if vcf_df.empty:
            logger.error("No variants matched the IDs from the ID list file")
            console.print("[red]Error:[/red] No variants matched the IDs from the ID list file")
            raise typer.Exit(code=1)

    logger.info(f"Loaded {len(vcf_df)} variants from VCF")

    logger.info("Loading reference genome...")
    # 使用 pyfaidx 加载参考基因组，as_raw=True 返回原始字节串提升性能
    fasta = Fasta(str(ref), as_raw=True)
    logger.info(f"Loaded reference with {len(fasta.keys())} chromosomes")

    out_list = []

    # 按染色体分组处理，避免重复访问
    for chrom in tqdm(vcf_df["chrom"].unique(), desc="Processing chromosomes"):
        if chrom not in fasta:
            logger.warning(f"Chromosome {chrom} not found in reference")
            continue

        logger.info(f"Processing chromosome: {chrom}")
        chrom_vcf_df = vcf_df[vcf_df["chrom"] == chrom]

        # 获取染色体序列对象
        chrom_seq = fasta[chrom]

        for row in tqdm(
            chrom_vcf_df.itertuples(),
            desc=f"Variants in {chrom}",
            total=len(chrom_vcf_df),
            leave=False,
        ):
            # 转换为0基坐标
            pos = row.pos - 1

            # 计算flanking区域的边界
            start = max(0, pos - flank_size)
            end = min(len(chrom_seq), pos + len(row.ref) + flank_size)

            # 使用pyfaidx的切片功能高效提取序列
            # as_raw=True 使得序列直接返回字符串，无需额外转换
            left_seq = chrom_seq[start:pos] if pos > start else ""
            right_seq = (
                chrom_seq[pos + len(row.ref) : end] if end > pos + len(row.ref) else ""
            )

            # 构建引物序列
            primer_seq = f"{left_seq}[{row.ref}/{row.alt}]{right_seq}"

            out_list.append({"name": row.effective_id, "sequence": primer_seq})

    logger.info(f"Generated {len(out_list)} primer sequences")

    # 输出结果
    out_df = pd.DataFrame(out_list)
    out_df.to_csv(out_file, sep="\t", index=False)
    logger.info(f"Results saved to {out_file}")


if __name__ == "__main__":
    typer.run(main)
