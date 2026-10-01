"""Tests for disease model detection"""
import pytest


def test_has_disease_model_patient():
    """Test detection of patient studies"""
    from scanner import _has_disease_model

    text = "We compared postmortem brain tissue from 20 schizophrenia patients."
    assert _has_disease_model(text) is True


def test_has_disease_model_genetic_mouse():
    """Test detection of genetic mouse models"""
    from scanner import _has_disease_model

    text = "APP/PS1 transgenic mice showed cognitive deficits."
    assert _has_disease_model(text) is True


def test_has_disease_model_chemically_induced():
    """Test detection of chemically-induced models"""
    from scanner import _has_disease_model

    text = "MPTP-treated mice were used as a Parkinson's disease model."
    assert _has_disease_model(text) is True


def test_has_disease_model_stress():
    """Test detection of stress models"""
    from scanner import _has_disease_model

    text = "Rats exposed to chronic mild stress showed depression-like behavior."
    assert _has_disease_model(text) is True


def test_has_disease_model_iPSC():
    """Test detection of iPSC-derived models"""
    from scanner import _has_disease_model

    text = "Patient-derived iPSC neurons were analyzed."
    assert _has_disease_model(text) is True


def test_has_disease_model_cell_lines_false():
    """Test that generic cell lines are NOT disease models"""
    from scanner import _has_disease_model

    text = "HEK293 cells were used for in vitro experiments."
    assert _has_disease_model(text) is False


def test_has_disease_model_no_model():
    """Test rejection when no disease model"""
    from scanner import _has_disease_model

    text = "Gene expression was analyzed using computational methods."
    assert _has_disease_model(text) is False


def test_has_disease_model_disorder_specific():
    """Test disorder-specific disease model patterns"""
    from scanner import _has_disease_model

    # Each disorder has specific patterns
    texts = [
        "schizophrenia patients showed cognitive deficits",
        "bipolar disorder subjects were compared to controls",
        "autism spectrum disorder children participated",
        "ADHD patients showed attention deficits",
        "major depression patients were enrolled",
        "anxiety disorder subjects completed the study",
        "PTSD patients were analyzed",
        "OCD patients participated",
        "Alzheimer's disease patients were studied",
        "Parkinson's disease patients were included",
    ]

    for text in texts:
        assert _has_disease_model(text) is True, f"Failed for: {text}"


def test_has_disease_model_knockout():
    """Test detection of knockout models"""
    from scanner import _has_disease_model

    text = "BDNF knockout mice showed memory deficits."
    assert _has_disease_model(text) is True


def test_has_disease_model_named_models():
    """Test detection of named genetic models"""
    from scanner import _has_disease_model

    # Named Alzheimer's models
    text = "5xFAD mice showed amyloid pathology."
    assert _has_disease_model(text) is True

    text = "APP-PS1 transgenic mice were used."
    assert _has_disease_model(text) is True
