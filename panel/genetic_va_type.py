import re


def is_missense_or_nonsynonymous(hgvs_p: str) -> str:
    if "*" in hgvs_p or "?" in hgvs_p:
        return "non-synonymous"
    protein_a, protein_b = re.split("[0-9]+", hgvs_p.split(".")[1])
    if protein_a == protein_b:
        return "synonymous"
    return "non-synonymous"
