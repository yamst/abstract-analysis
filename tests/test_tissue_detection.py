"""Tests for tissue detection functions"""
import pytest


def test_detect_tissue_brain():
    """Test detection of brain tissue"""
    from scanner import detect_tissue

    text = "Gene expression was measured in the prefrontal cortex."
    tissue = detect_tissue(text)

    assert tissue == "Brain"


def test_detect_tissue_gut():
    """Test detection of gut tissue"""
    from scanner import detect_tissue

    text = "We analyzed gene expression in the colon and intestine."
    tissue = detect_tissue(text)

    assert tissue == "Gut"


def test_detect_tissue_both():
    """Test detection of both tissues"""
    from scanner import detect_tissue

    text = "Brain and gut tissues were analyzed for gene expression."
    tissue = detect_tissue(text)

    assert tissue == "Brain & Gut"


def test_detect_tissue_other():
    """Test detection returns Other for no tissue keywords"""
    from scanner import detect_tissue

    text = "Gene expression was measured in the samples."
    tissue = detect_tissue(text)

    assert tissue == "Other"


def test_detect_tissue_for_gene_brain():
    """Test gene-specific tissue detection in brain"""
    from scanner import detect_tissue_for_gene

    text = "BDNF expression was decreased in the hippocampus."
    tissue = detect_tissue_for_gene("BDNF", text)

    assert tissue in ("Brain", "Brain & Gut")


def test_detect_tissue_for_gene_gut():
    """Test gene-specific tissue detection in gut"""
    from scanner import detect_tissue_for_gene

    text = "OCLN levels were altered in the intestinal epithelium."
    tissue = detect_tissue_for_gene("OCLN", text)

    assert tissue in ("Gut", "Brain & Gut")


def test_detect_tissue_brain_keywords():
    """Test various brain-related keywords"""
    from scanner import detect_tissue

    brain_keywords = [
        "prefrontal cortex",
        "hippocampus",
        "striatum",
        "cerebellum",
        "neuron",
        "glial",
        "synapse",
    ]

    for kw in brain_keywords:
        text = f"We measured gene expression in the {kw}."
        assert detect_tissue(text) == "Brain", f"Failed for keyword: {kw}"


def test_detect_tissue_gut_keywords():
    """Test various gut-related keywords"""
    from scanner import detect_tissue

    gut_keywords = [
        "colon",
        "intestine",
        "gut",
        "microbiome",
        "enteric",
        "fecal",
    ]

    for kw in gut_keywords:
        text = f"We analyzed {kw} samples."
        assert detect_tissue(text) == "Gut", f"Failed for keyword: {kw}"


def test_detect_tissue_case_insensitive():
    """Test that tissue detection is case-insensitive"""
    from scanner import detect_tissue

    assert detect_tissue("BRAIN tissue") == "Brain"
    assert detect_tissue("brain tissue") == "Brain"
    assert detect_tissue("Brain tissue") == "Brain"
