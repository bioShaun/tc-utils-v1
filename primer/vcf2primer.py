from pathlib import Path

import pandas as pd
import typer
from loguru import logger
from pyfaidx import Fasta
from tqdm import tqdm


def main(vcf_file: Path, ref: Path, out_file: Path, flank_size: int = 200) -> None:
    """Generate primer sequences from VCF and reference files."""
    logger.info("Reading VCF file...")
    vcf_df = pd.read_table(
        vcf_file,
        comment="#",
        usecols=[0, 1, 3, 4],
        names=["chrom", "pos", "ref", "alt"],
    )
    vcf_df["chrom"] = vcf_df["chrom"].astype("str")
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
        
        for row in tqdm(chrom_vcf_df.itertuples(), 
                       desc=f"Variants in {chrom}", 
                       total=len(chrom_vcf_df),
                       leave=False):
            # 转换为0基坐标
            pos = row.pos - 1
            
            # 计算flanking区域的边界
            start = max(0, pos - flank_size)
            end = min(len(chrom_seq), pos + len(row.ref) + flank_size)
            
            # 使用pyfaidx的切片功能高效提取序列
            # as_raw=True 使得序列直接返回字符串，无需额外转换
            left_seq = chrom_seq[start:pos] if pos > start else ""
            right_seq = chrom_seq[pos + len(row.ref):end] if end > pos + len(row.ref) else ""
            
            # 构建引物序列
            primer_seq = f"{left_seq}[{row.ref}/{row.alt}]{right_seq}"
            
            out_list.append({
                "name": f"{row.chrom}_{row.pos}", 
                "sequence": primer_seq
            })

    logger.info(f"Generated {len(out_list)} primer sequences")
    
    # 输出结果
    out_df = pd.DataFrame(out_list)
    out_df.to_csv(out_file, sep="\t", index=False)
    logger.info(f"Results saved to {out_file}")


if __name__ == "__main__":
    typer.run(main)    typer.run(main)