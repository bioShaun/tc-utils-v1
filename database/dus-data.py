def get_gt2(row):
    if row["genotype"] == "---":
        return "./."
    ref_count = row["genotype"].count(row["Ref Allele"])
    return {0: "1/1", 1: "0/1", 2: "0/0"}.get(ref_count)


def gt_seq(row):
    ref = row["ref"] + "/" + row["ref"]
    het = row["ref"] + "/" + row["alt"]
    alt = row["alt"] + "/" + row["alt"]
    return {"0/0": ref, "0/1": het, "1/1": alt, "./.": "---"}.get(row["final_gt"])
