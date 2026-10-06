# Known Issues

## Multi-Disorder Attribution Too Permissive

Abstracts are attributed to ALL mentioned disorders, not just studied ones.

**Example:** Schizophrenia postmortem study mentioned "cells lost earliest in Alzheimer's" as comparison → incorrectly attributed to Alzheimer's too.

**Fix:** Add `disorders_studied` field to LLM judgment (zero extra cost - same API call).
