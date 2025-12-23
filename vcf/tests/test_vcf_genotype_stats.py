def test_import():
    """Verify that the module can be imported."""
    try:
        from vcf import vcf_genotype_stats
    except ImportError:
        assert False, "Failed to import vcf.vcf_genotype_stats"
