"""Local-LLM qualification gate for scanner.py (hybrid mode).

Replaces the contextual qualification decision inside scanner.detect_gene_context
with a batched async call to a local llama.cpp / Ollama OpenAI-compatible server.
Gene extraction, per-gene tissue, and regulation stay as the cheap HGNC-validated
regex in scanner.py — the LLM only decides has_gene_context (qualifies: bool).

Run a server, e.g.:
    llama-server -m Qwen2.5-3B-Instruct-Q4_K_M.gguf -c 4096 -ngl 99 \
        -np 8 --cont-batching --host 127.0.0.1 --port 8080
or (smoke test only, slower, no continuous batching):
    OLLAMA_NUM_PARALLEL=4 ollama serve   # model: qwen2.5:3b-instruct-q4_K_M

scanner.py only calls qualify_batch_sync(items) -> {item_id: Judgment}.
"""

import os
import re
import json
import asyncio
import hashlib
import sqlite3
from dataclasses import dataclass

try:
    import httpx
except ImportError as e:
    raise ImportError(
        "llm_qualify requires httpx: pip install httpx"
    ) from e


# ── Config (env-overridable) ────────────────────────────────────────────────
MODEL_ID          = os.environ.get("LLM_MODEL", "qwen2.5-3b-instruct-q4_k_m")
PROMPT_VERSION    = os.environ.get("LLM_PROMPT_VERSION", "v1")   # bump to invalidate cache
LLM_BASE_URL      = os.environ.get("LLM_BASE_URL", "http://127.0.0.1:8080").rstrip("/")
LLM_CONCURRENCY   = int(os.environ.get("LLM_CONCURRENCY", "8"))
LLM_MAX_RETRIES   = int(os.environ.get("LLM_MAX_RETRIES", "2"))
LLM_TIMEOUT       = float(os.environ.get("LLM_TIMEOUT", "60"))
LLM_CACHE_PATH    = os.environ.get("LLM_CACHE_PATH", "llm_cache.sqlite")

_ENDPOINT = LLM_BASE_URL + "/v1/chat/completions"

# JSON schema the model must emit. We use response_format=json_object (broadly
# supported by llama.cpp and Ollama) and describe the schema in the prompt;
# _parse_judgment extracts defensively.
_RESPONSE_FORMAT = {"type": "json_object"}


@dataclass
class Judgment:
    qualifies: bool
    reason: str
    tissue_hint: str  # "brain" | "gut" | "both" | "none" ; advisory only


# ── Rubric (the exact inclusion criteria the LLM must judge) ─────────────────
SYSTEM_RUBRIC = """You are a strict biomedical-literature screener. Decide whether a PubMed abstract
QUALIFIES as a study we want, and return ONLY JSON.

A paper QUALIFIES (qualifies=true) ONLY IF ALL of these hold:
1. Real disease model: human patient / postmortem brain / clinical tissue, OR an animal
   disease model (genetic KO/KI/transgenic; or chemically-induced: MPTP, 6-OHDA, STZ,
   scopolamine, rotenone, cuprizone, kainic acid, poly-I:C; or stress models: chronic
   mild/unpredictable stress, social defeat, restraint, maternal separation), OR
   patient-derived iPSC neurons.
2. Differential gene expression, disease vs healthy/control: the paper reports gene
   expression that DIFFERS between a disease group and a healthy/control group
   (up/down-regulated, increased/decreased, differential expression in disease relative
   to control). "The gene is expressed in disease tissue" with no control comparison
   does NOT qualify.
3. Measured in brain and/or gut TISSUE: not serum/plasma/blood/CSF/urine/saliva only,
   not blood-derived cell lines (PBMC, lymphoblastoid/LCL) without tissue.
4. The change is attributable to the disorder itself, not to an intervention.

EXCLUDE (qualifies=false) if the expression change is caused by:
- a therapeutic drug / medication (antipsychotics, antidepressants, lithium, valproate,
  ketamine, stimulants, dexamethasone, etc.) administered to animals/patients;
- a non-drug intervention (DBS, TMS, tDCS, optogenetics, chemogenetics, exercise,
  environmental enrichment, diet/probiotic/prebiotic supplementation, fecal transplant,
  hyperbaric oxygen, surgery);
- expression attributed to a compound ("downregulated by X treatment", "X-mediated
  expression", "treatment altered expression");
- epigenetics / miRNA (methylation, histone marks, chromatin, miRNA/lncRNA/circRNA);
- behavioral-primary study (aim/primary result is behavioral/cognitive, expression only
  secondary);
- genetic-association / GWAS (polymorphism/SNP/genotype associated with disease risk,
  not expression levels);
- computational / bioinformatics reuse of public data (GEO, ArrayExpress, GTEx, PPI
  networks, in-silico) with no original measurement.

IMPORTANT distinctions (do not get these wrong):
- Chemically-induced DISEASE MODELS (MPTP, 6-OHDA, STZ, scopolamine, rotenone, cuprizone,
  kainic acid, poly-I:C) are VALID models, NOT exclusions. The gene change is attributable
  to the modeled disorder.
- Stress models (chronic mild stress, social defeat, restraint, maternal separation) are
  VALID depression/anxiety models, NOT intervention exclusions.
- "long-term depression" in a synaptic-plasticity context (with LTP, NMDA, AMPA, dendritic
  spines) is NOT a mood-disorder paper; judge the tissue/expression/intervention question
  on its merits.

You are given the abstract text plus a hint list of gene symbols already detected by a
regex/HGNC validator. Use them to locate the expression claim; DO NOT invent new gene
symbols and DO NOT reject just because a hint gene is absent.

Return JSON EXACTLY: {"qualifies": <bool>, "reason": "<one sentence citing which rule
fired>", "tissue_hint": "brain" | "gut" | "both" | "none"}
tissue_hint reflects where expression was measured in tissue (advisory; the authoritative
per-gene tissue is decided elsewhere). Output ONLY the JSON object."""

# ── Few-shot examples (hard cases) ───────────────────────────────────────────
_FEWSHOT = [
    # drug-induced (exclude) vs stress-induced (include)
    {
        "user": "Abstract: Haloperidol-treated mice showed decreased BDNF mRNA in the "
                "striatum versus vehicle controls. Detected genes: BDNF.",
        "assistant": json.dumps({"qualifies": False, "reason":
            "expression change caused by therapeutic drug haloperidol, not the disorder",
            "tissue_hint": "brain"}),
    },
    {
        "user": "Abstract: Rats exposed to chronic mild stress showed decreased hippocampal "
                "BDNF mRNA versus unstressed controls, modeling depression. Detected genes: BDNF.",
        "assistant": json.dumps({"qualifies": True, "reason":
            "chronic mild stress is a valid depression model; BDNF mRNA differs vs control "
            "in hippocampal tissue", "tissue_hint": "brain"}),
    },
    # behavioral-primary (exclude)
    {
        "user": "Abstract: We aimed to assess locomotor deficits in Shank3-KO mice; BDNF "
                "expression was also measured in the striatum. Detected genes: BDNF, SHANK3.",
        "assistant": json.dumps({"qualifies": False, "reason":
            "behavioral-primary study; gene expression is only secondary", "tissue_hint": "brain"}),
    },
    # GWAS (exclude) vs expression in a chemically-induced model (include)
    {
        "user": "Abstract: The DRD2 rs1800497 polymorphism was significantly associated with "
                "treatment response in schizophrenia patients. Detected genes: DRD2.",
        "assistant": json.dumps({"qualifies": False, "reason":
            "genetic-association/GWAS study of a polymorphism, not expression", "tissue_hint": "none"}),
    },
    {
        "user": "Abstract: MPTP-treated mice showed decreased DRD2 mRNA in the striatum "
                "versus saline controls, modeling Parkinson's disease. Detected genes: DRD2.",
        "assistant": json.dumps({"qualifies": True, "reason":
            "MPTP is a valid Parkinson model; DRD2 mRNA differs vs control in striatal tissue",
            "tissue_hint": "brain"}),
    },
    # fluid-only (exclude)
    {
        "user": "Abstract: Serum BDNF levels were decreased in depressed patients versus "
                "healthy controls; no tissue was measured. Detected genes: BDNF.",
        "assistant": json.dumps({"qualifies": False, "reason":
            "fluid-only measurement (serum); no brain/gut tissue", "tissue_hint": "none"}),
    },
]


def _build_messages(body, mentioned_genes):
    user = (
        f"Abstract:\n{body[:2500]}\n\n"
        f"Detected gene symbols (HGNC-validated, for reference): "
        f"{', '.join(sorted(mentioned_genes)) if mentioned_genes else '(none)'}\n\n"
        f"Return ONLY the JSON object."
    )
    msgs = [{"role": "system", "content": SYSTEM_RUBRIC}]
    for ex in _FEWSHOT:
        msgs.append({"role": "user", "content": ex["user"]})
        msgs.append({"role": "assistant", "content": ex["assistant"]})
    msgs.append({"role": "user", "content": user})
    return msgs


# ── SQLite cache (resume-safe, deterministic) ───────────────────────────────
def _cache_key(item_id):
    raw = f"{MODEL_ID}|{PROMPT_VERSION}|{item_id}"
    return hashlib.sha256(raw.encode()).hexdigest()


def _cache_connect():
    conn = sqlite3.connect(LLM_CACHE_PATH)
    conn.execute(
        """CREATE TABLE IF NOT EXISTS judgments (
            key            TEXT PRIMARY KEY,
            item_id        TEXT,
            qualifies      INTEGER,
            reason         TEXT,
            tissue_hint    TEXT,
            model_id       TEXT,
            prompt_version TEXT,
            ts             TEXT DEFAULT CURRENT_TIMESTAMP
        )"""
    )
    conn.execute("CREATE INDEX IF NOT EXISTS ix_judgments_item ON judgments(item_id)")
    conn.commit()
    return conn


def _cache_get_many(conn, keys):
    if not keys:
        return {}
    placeholders = ",".join("?" for _ in keys)
    cur = conn.execute(
        f"SELECT key, qualifies, reason, tissue_hint FROM judgments WHERE key IN ({placeholders})",
        keys,
    )
    out = {}
    for key, q, reason, th in cur.fetchall():
        out[key] = Judgment(bool(q), reason or "", th or "none")
    return out


def _cache_put_many(conn, rows):
    # rows: list[(key, item_id, Judgment)]
    if not rows:
        return
    conn.executemany(
        "INSERT OR REPLACE INTO judgments (key, item_id, qualifies, reason, tissue_hint, model_id, prompt_version) "
        "VALUES (?,?,?,?,?,?,?)",
        [(k, iid, int(j.qualifies), j.reason, j.tissue_hint, MODEL_ID, PROMPT_VERSION)
         for (k, iid, j) in rows],
    )
    conn.commit()


# ── Parsing ──────────────────────────────────────────────────────────────────
_JSON_OBJ_RE = re.compile(r"\{[\s\S]*\}")


def _parse_judgment(content):
    if not content:
        return Judgment(False, "empty_response", "none")
    try:
        obj = json.loads(content)
    except json.JSONDecodeError:
        m = _JSON_OBJ_RE.search(content)
        if not m:
            return Judgment(False, "json_parse_failed", "none")
        try:
            obj = json.loads(m.group(0))
        except json.JSONDecodeError:
            return Judgment(False, "json_parse_failed", "none")
    qualifies = bool(obj.get("qualifies", False))
    reason = str(obj.get("reason", "")).strip() or "no_reason"
    th = str(obj.get("tissue_hint", "none")).strip().lower()
    if th not in ("brain", "gut", "both", "none"):
        th = "none"
    return Judgment(qualifies, reason, th)


# ── Async client ─────────────────────────────────────────────────────────────
async def _qualify_one(client, sem, item_id, body, mentioned):
    payload = {
        "model": MODEL_ID,
        "messages": _build_messages(body, mentioned),
        "temperature": 0,
        "max_tokens": 96,
        "response_format": _RESPONSE_FORMAT,
        "stream": False,
    }
    last_err = None
    for attempt in range(LLM_MAX_RETRIES + 1):
        try:
            async with sem:
                resp = await client.post(_ENDPOINT, json=payload, timeout=LLM_TIMEOUT)
            if resp.status_code != 200:
                last_err = f"http_{resp.status_code}"
                await asyncio.sleep(0.5 * (attempt + 1))
                continue
            data = resp.json()
            content = data.get("choices", [{}])[0].get("message", {}).get("content", "")
            return item_id, _parse_judgment(content)
        except Exception as e:
            last_err = type(e).__name__
            await asyncio.sleep(0.5 * (attempt + 1))
    return item_id, Judgment(False, f"request_error:{last_err}", "none")


async def _qualify_batch(uncached, concurrency):
    sem = asyncio.Semaphore(concurrency)
    out = {}
    async with httpx.AsyncClient() as client:
        tasks = [_qualify_one(client, sem, iid, body, mentioned)
                 for (iid, body, mentioned) in uncached]
        for coro in asyncio.as_completed(tasks):
            iid, judgment = await coro
            out[iid] = judgment
    return out


# ── Public sync entry point (the only thing scanner.py calls) ────────────────
def qualify_batch_sync(items):
    """items: list[(item_id, body, mentioned_genes)] -> {item_id: Judgment}.

    item_id must be unique & stable (use the PMID when present, else a body hash).
    Serves cached judgments from SQLite; runs asyncio once for the uncached rest;
    persists new judgments in one transaction.
    """
    conn = _cache_connect()
    try:
        keys = {iid: _cache_key(iid) for (iid, _b, _m) in items}
        cached = _cache_get_many(conn, list(keys.values()))
        uncached = [(iid, body, mentioned)
                    for (iid, body, mentioned) in items
                    if keys[iid] not in cached]
        results = {iid: cached[keys[iid]] for iid in keys if keys[iid] in cached}

        if uncached:
            new = asyncio.run(_qualify_batch(uncached, LLM_CONCURRENCY))
            _cache_put_many(conn, [(keys[iid], iid, new[iid]) for iid in new])
            results.update(new)
    finally:
        conn.close()
    return results


if __name__ == "__main__":
    # Tiny self-test against a running server. Prints one judgment.
    sample = ("Rats exposed to chronic mild stress showed decreased hippocampal BDNF mRNA "
              "versus unstressed controls, modeling depression.")
    res = qualify_batch_sync([("selftest-1", sample, ["BDNF"])])
    j = res.get("selftest-1")
    print(j)
