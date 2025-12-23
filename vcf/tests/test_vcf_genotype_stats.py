import pytest
from vcf.vcf_genotype_stats import get_genotype_class

@pytest.mark.parametrize("a1, a2, expected", [
    (0, 0, "hom_ref"),
    (1, 1, "hom_alt"),
    (2, 2, "hom_alt"),
    (0, 1, "het"),
    (1, 0, "het"),
    (1, 2, "het"),
    (-1, -1, "missing"),
    (-1, 0, "missing"),
    (0, -1, "missing"),
    (1, -1, "missing"),
    (-1, 1, "missing"),
])
def test_get_genotype_class(a1, a2, expected):
    """Test the genotype classification logic."""
    assert get_genotype_class(a1, a2) == expected

def test_import():
    """Verify that the module can be imported."""
    try:
        from vcf import vcf_genotype_stats
    except ImportError:
        assert False, "Failed to import vcf.vcf_genotype_stats"