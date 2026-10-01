"""Tests for gene extraction from abstracts"""
import pytest


def test_detect_mentioned_genes_basic():
    """Test basic gene symbol detection"""
    from scanner import detect_mentioned_genes

    text = "BDNF expression was decreased in the hippocampus."
    genes = detect_mentioned_genes(text)

    assert "BDNF" in genes


def test_detect_mentioned_genes_multiple():
    """Test detection of multiple genes"""
    from scanner import detect_mentioned_genes

    text = "We measured BDNF, DRD2, and COMT mRNA levels in brain tissue."
    genes = detect_mentioned_genes(text)

    assert "BDNF" in genes
    assert "DRD2" in genes
    assert "COMT" in genes


def test_detect_mentioned_genes_alias():
    """Test detection via common aliases"""
    from scanner import detect_mentioned_genes

    # "tau" should resolve to MAPT
    text = "Tau protein levels were elevated in Alzheimer's patients."
    genes = detect_mentioned_genes(text)

    assert "MAPT" in genes


def test_detect_mentioned_genes_case_insensitive():
    """Test that gene detection is case-insensitive for aliases"""
    from scanner import detect_mentioned_genes

    # Mixed case should still work
    text = "BDNF and bdnf were measured."
    genes = detect_mentioned_genes(text)

    assert "BDNF" in genes


def test_detect_mentioned_genes_context_required():
    """Test that symbols need biological context"""
    from scanner import detect_mentioned_genes

    # "CAT" alone might not be detected without context
    text = "The CAT scan showed abnormalities."
    genes = detect_mentioned_genes(text)

    # CAT without biological context should not be detected
    assert "CAT" not in genes


def test_detect_mentioned_genes_excludes_non_genes():
    """Test that non-gene acronyms are excluded"""
    from scanner import detect_mentioned_genes

    # "APP" is a gene but "AP" alone is not typically
    text = "AP was used to analyze the data."
    genes = detect_mentioned_genes(text)

    # Should not detect AP as a gene
    assert "AP" not in genes


def test_detect_mentioned_genes_in_sentence():
    """Test gene detection in a full sentence context"""
    from scanner import detect_mentioned_genes

    text = """
    Schizophrenia patients showed altered gene expression.
    BDNF mRNA was significantly decreased in prefrontal cortex
    compared to healthy controls. DRD2 levels were unchanged.
    """
    genes = detect_mentioned_genes(text)

    assert "BDNF" in genes
    assert "DRD2" in genes
