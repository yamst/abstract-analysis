"""Tests for regex pre-filter and exclusion predicates"""
import pytest


def test_abstract_uses_therapeutic_drug():
    """Test that therapeutic drug studies are excluded"""
    from scanner import _abstract_uses_therapeutic_drug

    # Should exclude: drug administration
    text = "Mice were treated with haloperidol for 14 days."
    assert _abstract_uses_therapeutic_drug(text) is True

    # Should exclude: medication mentioned
    text = "Patients receiving fluoxetine showed improvement."
    assert _abstract_uses_therapeutic_drug(text) is True

    # Should NOT exclude: no drug
    text = "Gene expression was measured in postmortem brain tissue."
    assert _abstract_uses_therapeutic_drug(text) is False


def test_abstract_uses_epigenetics():
    """Test that epigenetics studies are excluded"""
    from scanner import _abstract_uses_epigenetics

    # Should exclude: miRNA
    text = "We measured miR-132 expression in the hippocampus."
    assert _abstract_uses_epigenetics(text) is True

    # Should exclude: DNA methylation
    text = "DNA methylation patterns were altered in patients."
    assert _abstract_uses_epigenetics(text) is True

    # Should NOT exclude: regular gene expression
    text = "BDNF mRNA levels were decreased in the prefrontal cortex."
    assert _abstract_uses_epigenetics(text) is False


def test_is_genetic_association_study():
    """Test that GWAS/genetic association studies are excluded"""
    from scanner import _is_genetic_association_study

    # Should exclude: GWAS
    text = "A genome-wide association study identified risk variants."
    assert _is_genetic_association_study(text) is True

    # Should exclude: SNP association
    text = "The rs1234 polymorphism was associated with disease risk."
    assert _is_genetic_association_study(text) is True

    # Should NOT exclude: expression study
    text = "Gene expression levels were measured using RNA-seq."
    assert _is_genetic_association_study(text) is False


def test_regex_prefilter_passes():
    """Test cases that should pass the prefilter"""
    from scanner import _regex_prefilter

    # Valid patient study with expression measurement
    text = """
    BACKGROUND: Schizophrenia is associated with cognitive dysfunction.
    METHODS: We compared postmortem brain tissue from 20 schizophrenia patients
    and 20 healthy controls.
    RESULTS: BDNF mRNA expression was significantly decreased in the prefrontal
    cortex of patients compared to controls (p<0.01).
    """
    assert _regex_prefilter(text) is True


def test_regex_prefilter_excludes_drugs():
    """Test that drug studies fail prefilter"""
    from scanner import _regex_prefilter

    text = """
    BACKGROUND: Schizophrenia is treated with antipsychotics.
    METHODS: Mice were treated with clozapine for 21 days.
    RESULTS: Clozapine treatment altered BDNF expression in the striatum.
    """
    assert _regex_prefilter(text) is False


def test_regex_prefilter_excludes_gwas():
    """Test that GWAS studies fail prefilter"""
    from scanner import _regex_prefilter

    text = """
    BACKGROUND: Genetic factors contribute to Alzheimer's disease risk.
    METHODS: We performed a genome-wide association study in 5000 patients.
    RESULTS: The APOE rs429358 variant was strongly associated with AD risk.
    """
    assert _regex_prefilter(text) is False
