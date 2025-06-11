#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from pathlib import Path
from typing import List

import pandas as pd
import typer
from loguru import logger

app = typer.Typer()


def parse_fastq_filename(sample_path: Path):
    filename = sample_path.name
    fastqs = list(sample_path.glob("*.fastq.gz"))

    lib_id = filename.split("-")[-1]
    lib_info = []

    for fq in fastqs:
        if fq.name.endswith("combined_R1.fastq.gz"):
            read_type = "R1"
        elif fq.name.endswith("combined_R2.fastq.gz"):
            read_type = "R2"
        else:
            raise Exception("无法识别的fastq文件名")
        lib_info.append(
            {"libid": lib_id, "read_type": read_type, "path": fq.absolute()}
        )
    return lib_info


def build_libid_fastq_map(fastq_path: Path) -> pd.DataFrame:
    libid_map = []
    for path in fastq_path.glob("Sample*"):
        try:
            info = parse_fastq_filename(path)
            libid_map.extend(info)
        except Exception:
            continue
    return pd.DataFrame(libid_map)


def read_or_build_config(fq_line_dir: Path) -> pd.DataFrame:
    config_file = fq_line_dir / "libid_fastq_config.tsv"
    if config_file.exists():
        typer.echo(f"使用已有配置文件: {config_file}")
        return pd.read_table(config_file)

    libid_map = build_libid_fastq_map(fq_line_dir)
    libid_map["dir_name"] = fq_line_dir.name
    libid_map.to_csv(config_file, sep="\t", index=False)
    return libid_map


def load_config(base_dir: Path) -> pd.DataFrame:
    libid_map_list = []
    for each_path in base_dir.glob("*/tcwl-*"):
        if each_path.parent.name.startswith("20"):
            logger.info(f"获取libid-fastq配置：{each_path.name}")
            libid_map = read_or_build_config(each_path)
            libid_map_list.append(libid_map)
    all_libid_map = pd.concat(libid_map_list)
    return all_libid_map


def group_fastq_by_sample(libid_map, libid_to_sample):
    sample_files = defaultdict(lambda: defaultdict(list))
    for libid, reads in libid_map.items():
        sample_id = libid_to_sample.get(libid)
        if sample_id:
            for read_type, paths in reads.items():
                sample_files[sample_id][read_type].extend(paths)
    return sample_files


def merge_or_link_sh(fq_list: List[str], name: str):
    if len(fq_list) == 1:
        return f"ln -s {fq_list[0]} {name}"
    return f"cat {' '.join(fq_list)} > {name}"


def write_nextflow_input(fq_df: pd.DataFrame, output_dir: Path):
    for (sample_id, read_type), sample_df in fq_df.groupby(["sample_id", "read_type"]):
        out_fq = f"{sample_id}.{read_type}.fq.gz"
        fq_list = sorted(sample_df["path"].to_list())
        cmd_file = output_dir / f"mergeFastq-{sample_id}-{read_type}.sh"
        cmd = merge_or_link_sh(fq_list, out_fq)
        with open(cmd_file, "w") as f:
            f.write(f"{cmd}\n")


@app.command()
def main(
    sample_info: Path = typer.Option(
        ..., help="样品信息TSV文件，必须包含libid、sample_id、batch_dir列"
    ),
    base_dir: Path = typer.Option(..., help="包含所有FASTQ路径的文本文件"),
    data_cuoff: float = typer.Option(0.01, help="数据量阈值"),
    output_dir: Path = typer.Option(
        "nextflow_dir", help="Nextflow输入bash文件保存目录"
    ),
):
    sample_df = pd.read_table(
        sample_info,
        header=None,
        names=["libid", "sample_id", "data_size", "dir_name"],
        usecols=[0, 1, 2, 3],
    )

    passed_sample_df = sample_df[sample_df["data_size"] >= data_cuoff]

    libid_map = load_config(base_dir)
    add_fq_df = passed_sample_df.merge(libid_map, how="left")
    # add_fq_df.to_csv(output, sep="\t", index=False)
    output_dir.mkdir(exist_ok=True, parents=True)
    write_nextflow_input(add_fq_df, output_dir)
    typer.echo(f"已写入Nextflow输入文件: {output_dir}")


if __name__ == "__main__":
    app()
