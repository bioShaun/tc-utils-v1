import re
from pathlib import Path

import typer


def main(gtf: Path) -> None:
    out_gtf = gtf.with_suffix(".addGene.gtf")
    tr_dict = {}
    with gtf.open("r") as in_gtf:
        for line in in_gtf:
            line_inf = line.strip().split("\t")
            out_inf = line
            # tr_id = line_inf[-1].split()[1]
            tr_id = re.search(r'transcript_id "(\S+)"', line_inf[-1]).groups()[0]
            strand = line_inf[6]
            if strand not in ["+", "-"]:
                line_inf[6] = "+"
            if tr_id not in tr_dict:
                tr_dict[tr_id] = {}
                tr_dict[tr_id]["exon"] = []
                tr_dict[tr_id]["CDS"] = []
                tr_dict[tr_id]["transcript"] = []
            new_out_inf = "\t".join(line_inf)
            if "gene_id" not in line:
                gene_id = line.strip().split()[-1]
                new_out_inf = f"{line.strip()} gene_id {gene_id}\n"
            tr_dict[tr_id][line_inf[2]].append(new_out_inf + "\n")
    with out_gtf.open("w") as out:
        for transcript in tr_dict:
            for feature in ["transcript", "exon", "CDS"]:
                feature_lines = tr_dict[transcript][feature]
                print(feature_lines)
                if len(feature_lines) == 0:
                    print(transcript, feature)
                    if feature == "exon":
                        feature_lines = tr_dict[transcript]["CDS"]
                    else:
                        feature_lines = tr_dict[transcript]["exon"]
                for line in feature_lines:
                    line_list = line.split("\t")
                    line_list[2] = feature
                    out.write("\t".join(line_list))


if __name__ == "__main__":
    typer.run(main)
