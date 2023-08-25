from pathlib import Path
from typing import Dict
import typer


def gff_tag_dict(gff_info: str) -> Dict[str, str]:
    gff_dict = {}
    gff_list = gff_info.split(";")
    for item in gff_list:
        key, value = item.split("=")
        gff_dict[key] = value
    return gff_dict


def main(gff: Path) -> None:
    with gff.open() as gff_inf:
        for eachline in gff_inf:
            eachline_list = eachline.strip().split("\t")
            if eachline_list[2] == "mRNA":
                eachline_dict = gff_tag_dict(eachline_list[-1])
                each_name = eachline_dict.get("Name", "")
                each_id = eachline_dict.get("Parent", "")
                print(
                    f"{eachline_list[0]}\t{eachline_list[3]}\t{eachline_list[4]}\t{each_id}\t{each_name}"
                )


if __name__ == "__main__":
    typer.run(main)
