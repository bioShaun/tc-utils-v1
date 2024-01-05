import re
from pathlib import Path
from typing import List, Optional

import delegator
import pandas as pd
import typer
from tqdm import tqdm

DATA_COLS = ["accession", "name", "dateSize", "deliverDir", "storeDir", "libId"]
LAUNCH_SH_PATH = "/public/home/zxchen/scripts/system/omsrun_one.sh"
LAUNCH_SH_PATH_NEW = "/public/home/zxchen/scripts/system/omsrun_one_new.sh"


def launch_cmd(cmd_path: Optional[Path], new_queue: bool = False) -> None:
    if cmd_path:
        if new_queue:
            delegator.run(f"{LAUNCH_SH_PATH_NEW} {cmd_path}")
        else:
            delegator.run(f"{LAUNCH_SH_PATH} {cmd_path}")


def write_cmd(
    fq_list: List[str], name: str, number: int, out_dir: Path, mode="link"
) -> Optional[Path]:
    fq_str_list = [str(each) for each in fq_list]
    cmd = ""
    if len(fq_str_list) == 1:
        fq_out = out_dir / f"{name}.R{number}.fq.gz"
        if not fq_out.is_file():
            if mode == "link":
                cmd = f'ln -s {" ".join(fq_str_list)} {out_dir}/{name}.R{number}.fq.gz'
                delegator.run(cmd)
                return None
            elif mode == "cp":
                cmd = f'cp {" ".join(fq_str_list)} {out_dir}/{name}.R{number}.fq.gz'
            else:
                raise ValueError(f"mode [{mode}] not allowed")
        else:
            return None
    else:
        cmd = f'cat {" ".join(fq_str_list)} > {out_dir}/{name}.R{number}.fq.gz'
    script_dir = out_dir / "scripts"
    script_dir.mkdir(exist_ok=True, parents=True)
    script_path = script_dir / f"mergeFastq-{name}-{number}.sh"
    with open(script_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"{cmd}\n")
    return script_path


def file_patterns(deliverDir, storeDir, libname):
    patterns = [
        f"*/{deliverDir}/*{storeDir}*/**/*{libname}*/*_{libname}_*.",
        f"*/{deliverDir}/*{storeDir}*/**/*_{libname}_*.",
        f"*/{deliverDir}/*{storeDir}*/**/*{libname}*/*{libname}_*.",
        f"*/{deliverDir}/*{storeDir}*/**/*{libname}_*.",
        f"*/{deliverDir}/*{storeDir}*/**/*{libname}-*.",
    ]
    fq_patterns = []
    for pattern in patterns:
        for fq_label in ["fq", "fastq"]:
            fq_patterns.append(f"{pattern}{fq_label}.gz")
    return fq_patterns


def main(
    data_path: Path,
    sample_file: Path,
    out_dir: Path,
    generate: bool = typer.Option(False, "--generate"),
    run: bool = typer.Option(False, "--run"),
    mode: str = typer.Option("link", "--mode"),
    fq_label: str = "fq",
    new_queue: bool = False,
) -> None:
    sample_df = pd.read_csv(sample_file, sep="\t", header=None, names=DATA_COLS)
    dup_df = sample_df[
        sample_df.duplicated(subset=["libId", "deliverDir", "storeDir"])
    ].copy()
    if not dup_df.empty:
        print(dup_df)
        raise ValueError("Duplicated library.")
    sample_df["name"] = sample_df["name"].map(
        lambda x: re.sub("[^a-zA-Z0-9_-]", "_", str(x))
    )
    for name, df in tqdm(sample_df.groupby("name")):
        fq_files = []
        for index, row in df.iterrows():
            libId = row["libId"]
            deliverDir = row["deliverDir"]
            storeDir = row["storeDir"]
            # fq_dir = data_path / deliverDir / storeDir
            if libId.startswith("CWHPE"):
                libname = "-".join(libId.split("-")[1:])
            else:
                libname = libId
            for pattern in file_patterns(deliverDir, storeDir, libname):
                # print(pattern)
                lib_fq_files = list(data_path.glob(pattern))
                if lib_fq_files:
                    fq_files.extend(lib_fq_files)
                    break
            if len(lib_fq_files) == 0:
                print(f"No fq files found for {libId} {deliverDir}-{storeDir}")

        fq_files.sort()
        fq1 = [each for each in fq_files if re.search("1.(fq|fastq).gz$", each.name)]
        fq2 = [each for each in fq_files if re.search("2.(fq|fastq).gz$", each.name)]
        fq_script1 = write_cmd(fq1, name, 1, out_dir, mode)
        fq_script2 = write_cmd(fq2, name, 2, out_dir, mode)
        if run:
            if not (fq_script1 is None):
                launch_cmd(fq_script1, new_queue)
            if not (fq_script2 is None):
                launch_cmd(fq_script2, new_queue)


if __name__ == "__main__":
    typer.run(main)
