import pandas as pd

from chip.gt_to_pheno_update import format_output


def test_format_output_adds_lc_suffix_only() -> None:
    input_df = pd.DataFrame(
        {
            "location_status": ["HC", "LC"],
            "CHROM": ["1A", "1B"],
            "POS": [100, 200],
            "category": ["产量", "产量"],
            "trait": ["千粒重", "千粒重"],
            "gene": ["GeneHC", "GeneLC"],
            "sample1": ["高千粒重", "低千粒重"],
        }
    )

    result_df = format_output(input_df)

    assert ("千粒重", "GeneHC") in result_df.columns
    assert ("千粒重", "GeneLC(LC)") in result_df.columns
    assert ("千粒重", "GeneLC") not in result_df.columns
    assert result_df.loc["sample1", ("千粒重", "GeneHC")] == "高千粒重"
    assert result_df.loc["sample1", ("千粒重", "GeneLC(LC)")] == "低千粒重"
