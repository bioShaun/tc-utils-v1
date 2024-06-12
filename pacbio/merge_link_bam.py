from pathlib import Path

import delegator
import typer


def main(pacbio_dir: Path, out_dir: Path) -> None:
    pacbio_abs_dir = pacbio_dir.resolve()
    for each in pacbio_abs_dir.iterdir():
        bams = list(each.glob("**/*.bam"))
        if len(bams) == 1:
            cmd = f"ln -s {bams[0]} {out_dir}/{each.name}.bam"
            delegator.run(cmd)
        else:
            bams_str = " ".join([str(each) for each in bams])
            cmd = f"pbmerge {bams_str} > {out_dir}/{each.name}.bam"
            print(cmd)


if __name__ == "__main__":
    typer.run(main)
