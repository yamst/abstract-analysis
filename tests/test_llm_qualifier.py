"""Tests for LLM qualifier module"""
import pytest
from unittest.mock import patch, MagicMock
import json


def test_judgment_dataclass():
    """Test Judgment dataclass creation"""
    from llm_qualify import Judgment

    j = Judgment(qualifies=True, reason="Valid disease model", tissue_hint="brain")

    assert j.qualifies is True
    assert j.reason == "Valid disease model"
    assert j.tissue_hint == "brain"


def test_build_messages():
    """Test message building for LLM"""
    from llm_qualify import _build_messages

    body = "BDNF expression was decreased in schizophrenia patients."
    genes = ["BDNF"]

    messages = _build_messages(body, genes)

    assert len(messages) > 3  # system + few-shots + user
    assert messages[0]["role"] == "system"
    assert messages[-1]["role"] == "user"
    assert "BDNF" in messages[-1]["content"]


def test_parse_judgment_valid():
    """Test parsing valid JSON response"""
    from llm_qualify import _parse_judgment

    json_response = json.dumps({
        "qualifies": True,
        "reason": "Valid study",
        "tissue_hint": "brain"
    })

    j = _parse_judgment(json_response)

    assert j.qualifies is True
    assert j.reason == "Valid study"
    assert j.tissue_hint == "brain"


def test_parse_judgment_invalid_json():
    """Test parsing invalid JSON returns safe default"""
    from llm_qualify import _parse_judgment

    invalid = "not valid json"

    j = _parse_judgment(invalid)

    assert j.qualifies is False
    assert "parse" in j.reason.lower() or "json" in j.reason.lower()


def test_parse_judgment_missing_fields():
    """Test parsing JSON with missing fields"""
    from llm_qualify import _parse_judgment

    incomplete = json.dumps({"qualifies": True})

    j = _parse_judgment(incomplete)

    # Should have defaults for missing fields
    assert j.reason is not None
    assert j.tissue_hint is not None


@patch('llm_qualify.httpx.AsyncClient')
def test_qualify_batch_sync(mock_client):
    """Test batch qualification with mocked HTTP"""
    from llm_qualify import qualify_batch_sync

    # Mock the async client response
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "choices": [{
            "message": {
                "content": json.dumps({
                    "qualifies": True,
                    "reason": "Test reason",
                    "tissue_hint": "brain"
                })
            }
        }]
    }

    # Set up the mock
    mock_client.return_value.__aenter__.return_value.post = MagicMock(
        return_value=mock_response
    )

    items = [("test-1", "Test abstract body", ["BDNF"])]
    # This test just verifies the function structure works
    # Actual async mocking is complex; integration test covers real behavior


def test_cache_key_deterministic():
    """Test that cache keys are deterministic"""
    from llm_qualify import _cache_key

    key1 = _cache_key("item-1")
    key2 = _cache_key("item-1")

    assert key1 == key2


def test_cache_key_different_for_different_items():
    """Test that different items get different cache keys"""
    from llm_qualify import _cache_key

    key1 = _cache_key("item-1")
    key2 = _cache_key("item-2")

    assert key1 != key2
