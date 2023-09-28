import re
from pathlib import Path
from typing import Dict, List

import typer
from attrs import define, field

app = typer.Typer()


def get_tag_value_from_gtf(tag: str, info: str) -> str:
    pattern = f'{tag} "(.*?)";'
    if re.findall(pattern, info):
        return re.findall(pattern, info)[0]
    else:
        return ""


@define
class Gene:
    gene_id: str
    name: str
    out_id: str
    refseq_id: str = field(default=None)
    position: str = field(init=False)
    description_set: set[str] = field(init=False)

    def __attrs_post_init__(self):
        self.description_set = set()

    def add_info(self, gtf_line: str) -> None:
        gene_inf = gtf_line.strip().split("\t")
        gene_attrs = gene_inf[8]
        if gene_inf[2] == "gene":
            self.add_position(gene_inf)
        else:
            description = get_tag_value_from_gtf("product", gene_attrs)
            self.add_description(description)
        refseq_id = get_tag_value_from_gtf("db_xref", gene_attrs)
        if refseq_id:
            self.refseq_id = refseq_id.split(":")[1]

    def add_position(self, gene_inf: List[str]) -> None:
        self.position = f"{gene_inf[0]}:{gene_inf[3]}-{gene_inf[4]}|{gene_inf[6]}"

    def add_description(self, description: str) -> None:
        if description:
            self.description_set.add(description)

    @property
    def gene_out_id(self) -> str:
        if self.out_id == "refseq_id":
            if self.refseq_id:
                return self.refseq_id
        return self.gene_id

    @property
    def description(self) -> str:
        if self.description_set:
            return " | ".join(self.description_set)
        else:
            return "--"


@define
class Transcript:
    tr_id: str
    gene_id: str
    protein_id: str = field(default=None)

    def add_info(self, gtf_line: str) -> None:
        gene_inf = gtf_line.strip().split("\t")
        gene_attrs = gene_inf[8]
        protein_id = get_tag_value_from_gtf("protein_id", gene_attrs)
        if protein_id:
            self.protein_id = protein_id


@app.command()
def gene_annotation(
    gtf: Path, annotation: Path, gene_id_att: str = "refseq_id"
) -> None:
    gene_dict: dict[str, Gene] = {}
    with open(gtf) as gtf_inf:
        for eachline in gtf_inf:
            if eachline.startswith("#"):
                continue
            eachline_inf = eachline.strip().split("\t")
            eachline_attrs = eachline_inf[8]
            gene_id = get_tag_value_from_gtf("gene_id", eachline_attrs)
            gene_name = get_tag_value_from_gtf("gene", eachline_attrs)
            if gene_id not in gene_dict:
                gene_dict[gene_id] = Gene(
                    gene_id=gene_id, name=gene_name, out_id=gene_id_att
                )
            gene_dict[gene_id].add_info(eachline)

    with open(annotation, "w") as annotation_inf:
        annotation_inf.write("gene_id\tLocus\tgene_name\tgene_description\n")
        for each_gene in gene_dict.values():
            try:
                annotation_inf.write(
                    f"{each_gene.gene_out_id}\t{each_gene.position}\t{each_gene.name}\t{each_gene.description}\n"
                )
            except AttributeError:
                print(each_gene)
                raise AttributeError()


@app.command()
def transcript_annotation(gtf: Path, annotation: Path) -> None:
    tr_dict: dict[str, Transcript] = {}
    gene_dict: dict[str, str] = {}
    with open(gtf) as gtf_inf:
        for eachline in gtf_inf:
            if eachline.startswith("#"):
                continue
            eachline_inf = eachline.strip().split("\t")
            eachline_attrs = eachline_inf[8]
            gene_id = get_tag_value_from_gtf("gene_id", eachline_attrs)
            if eachline_inf[2] == "gene":
                refseq_id = get_tag_value_from_gtf("db_xref", eachline_attrs)
                gene_dict[gene_id] = refseq_id.split(":")[1]
                continue

            tr_id = get_tag_value_from_gtf("transcript_id", eachline_attrs)

            if tr_id not in tr_dict:
                tr_dict[tr_id] = Transcript(tr_id=tr_id, gene_id=gene_id)
            tr_dict[tr_id].add_info(eachline)

    with open(annotation, "w") as annotation_inf:
        for each_tr in tr_dict.values():
            try:
                refseq_id = gene_dict[each_tr.gene_id]
                annotation_inf.write(
                    f"{refseq_id}\t{each_tr.tr_id}\t{each_tr.protein_id}\n"
                )
            except:
                print(each_tr)
                raise AttributeError()


if __name__ == "__main__":
    app()
