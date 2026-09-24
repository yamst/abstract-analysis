import re
import socket
socket.setdefaulttimeout(60)   # abort any PubMed request that hangs > 60 s
from Bio import Entrez
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from collections import defaultdict
import time
import datetime
import os
import hashlib
import ssl
import urllib.request


_HGNC_FULL_FILE = os.path.join(os.path.dirname(os.path.abspath(globals().get('__file__', 'scanner.py'))), "hgnc_full.txt")
_HGNC_FULL_URL  = ("https://www.genenames.org/cgi-bin/download/custom?"
                   "col=gd_app_sym&col=gd_prev_sym&col=gd_aliases"
                   "&status=Approved&hgnc_dbtag=on"
                   "&order_by=gd_app_sym_sort&format=text&submit=submit")

def _load_hgnc():
    if not os.path.exists(_HGNC_FULL_FILE):
        print("Downloading HGNC gene symbol + alias list…")
        urllib.request.urlretrieve(_HGNC_FULL_URL, _HGNC_FULL_FILE)
    symbols = set()
    alias_to_sym = {}   # uppercase alias → approved symbol
    with open(_HGNC_FULL_FILE) as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            approved = parts[0].strip()
            if not approved:
                continue
            symbols.add(approved.upper())
            prev    = parts[1].strip() if len(parts) > 1 else ""
            aliases = parts[2].strip() if len(parts) > 2 else ""
            for raw in [prev, aliases]:
                for alias in raw.split(","):
                    alias = alias.strip()
                    if alias and alias != approved:
                        alias_to_sym[alias.upper()] = approved
    print(f"HGNC loaded: {len(symbols):,} approved symbols, {len(alias_to_sym):,} aliases")
    return symbols, alias_to_sym

HGNC_SYMBOLS, _HGNC_ALIAS_TO_SYM = _load_hgnc()



ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "dasfour@outlook.com"


# ── LLM qualification gate (hybrid mode) ────────────────────────────────────
# When USE_LLM_QUALIFIER=1, the contextual qualification decision currently inside
# detect_gene_context (the ~11 _abstract_* exclusion predicates + the disease-model
# and expression-measurement positive requirements) is replaced by a batched local
# LLM call (see llm_qualify.py). Gene extraction, per-gene tissue attribution, and
# regulation direction stay as the cheap HGNC-validated regex below; the LLM only
# decides has_gene_context. The regex gate in detect_gene_context is retained as the
# fallback path (server down) and the eval baseline (eval_qualify.py).
USE_LLM_QUALIFIER = os.environ.get("USE_LLM_QUALIFIER", "0") == "1"
_llm_qualify_fn = None
if USE_LLM_QUALIFIER:
    try:
        from llm_qualify import qualify_batch_sync as _llm_qualify_fn
    except Exception as _e:
        print(f"[warn] USE_LLM_QUALIFIER=1 but llm_qualify unavailable ({_e}); "
              f"falling back to the regex gate.")
        USE_LLM_QUALIFIER = False


DISORDERS = {
    "Schizophrenia":    (
        "Schizophrenia[MeSH Terms] OR schizophrenia[Title/Abstract] OR "
        "schizophrenic[Title/Abstract] OR schizoaffective[Title/Abstract]"
    ),
    "Bipolar Disorder": (
        "Bipolar Disorder[MeSH Terms] OR bipolar disorder[Title/Abstract] OR "
        "bipolar depression[Title/Abstract] OR manic depression[Title/Abstract] OR "
        "manic episode[Title/Abstract]"
    ),
    "Autism":           (
        "Autistic Disorder[MeSH Terms] OR Autism Spectrum Disorder[MeSH Terms] OR "
        "autism[Title/Abstract] OR autistic[Title/Abstract] OR "
        "asperger[Title/Abstract]"
    ),
    "ADHD":             (
        "Attention Deficit Disorder with Hyperactivity[MeSH Terms] OR "
        "attention deficit hyperactivity disorder[Title/Abstract] OR "
        "ADHD[Title/Abstract] OR attention deficit disorder[Title/Abstract]"
    ),
    "Major Depression": (
        "Depressive Disorder, Major[MeSH Terms] OR Depressive Disorder[MeSH Terms] OR "
        "major depressive disorder[Title/Abstract] OR major depression[Title/Abstract] OR "
        "unipolar depression[Title/Abstract]"
    ),
    "Anxiety":          (
        "Anxiety Disorders[MeSH Terms] OR Anxiety[MeSH Terms] OR "
        "anxiety disorder[Title/Abstract] OR anxiety[Title/Abstract] OR "
        "generalized anxiety[Title/Abstract] OR social anxiety[Title/Abstract]"
    ),
    "PTSD":             (
        "Stress Disorders, Post-Traumatic[MeSH Terms] OR "
        "post-traumatic stress disorder[Title/Abstract] OR PTSD[Title/Abstract] OR "
        "posttraumatic stress[Title/Abstract] OR post traumatic stress[Title/Abstract]"
    ),
    "OCD":              (
        "Obsessive-Compulsive Disorder[MeSH Terms] OR "
        "obsessive compulsive disorder[Title/Abstract] OR OCD[Title/Abstract] OR "
        "obsessive-compulsive[Title/Abstract]"
    ),
    "Alzheimer's":      (
        "Alzheimer Disease[MeSH Terms] OR alzheimer disease[Title/Abstract] OR "
        "alzheimer's disease[Title/Abstract] OR alzheimer[Title/Abstract]"
    ),
    "Parkinson's":      (
        "Parkinson Disease[MeSH Terms] OR parkinson disease[Title/Abstract] OR "
        "parkinson's disease[Title/Abstract] OR parkinsonism[Title/Abstract]"
    ),
}

RODENT_SPECIES = {
    "Mice":            ["mouse", "mice", "murine", "mus musculus"],
    "Rats":            ["rat", "rats", "rattus norvegicus", "rattus"],
    "Gerbils":         ["gerbil", "gerbils", "meriones"],
    "Guinea Pigs":     ["guinea pig", "guinea pigs", "cavia porcellus"],
    "Hamsters":        ["hamster", "hamsters", "cricetulus", "mesocricetus"],
    "Voles":           ["vole", "voles", "microtus"],
    "Squirrels":       ["squirrel", "squirrels", "sciurus"],
    "Prairie Dogs":    ["prairie dog", "prairie dogs", "cynomys"],
    "Chinchillas":     ["chinchilla", "chinchillas"],
    "Degus":           ["degu", "degus", "octodon"],
    "Naked Mole-Rats": ["naked mole-rat", "naked mole rat", "heterocephalus glaber"],
}

SPECIES_KEYWORDS = {
    "Humans":     ["human subjects", "human participants", "healthy volunteers",
                   "clinical cohort", "human patients", "human", "humans",
                   "IPSC", "induced pluripotent stem cell", "iPSC",
                   "brain organoid", "cerebral organoid", "clinical trial",
                   "postmortem tissue", "post mortem tissue"],
    "Drosophila": ["drosophila", "fruit fly", "drosophila melanogaster",],
    "C. elegans": ["c. elegans", "caenorhabditis elegans"],
    "Zebrafish":  ["zebrafish", "danio rerio"],
    "Monkeys":    ["monkey", "monkeys", "macaque", "primate", "non-human primate"],
    "Dogs":       ["dog", "dogs", "canine"],
    "Pigs":       ["pig", "pigs", "porcine"],
    "Yeast":      ["yeast", "saccharomyces"],
    "Frogs":      ["frog", "frogs", "xenopus"],
}

BRAIN_KEYWORDS = [
    "brain", "cortex", "hippocampus", "hippocampal", "striatum", "striatal",
    "prefrontal", "amygdala", "cerebellum", "cerebellar", "neuron", "neurons",
    "neuronal", "synapse", "synaptic", "dopamine", "serotonin", "gaba",
    "glutamate", "neural", "neurological", "cerebral", "basal ganglia",
    "substantia nigra", "thalamus", "hypothalamus", "brainstem", "spinal cord",
    "astrocyte", "microglia", "oligodendrocyte", "dendrite", "axon",
]

GUT_KEYWORDS = [
    "gut", "intestine", "intestinal", "colon", "colonic", "gastrointestinal",
    "microbiome", "microbiota", "gut-brain", "gut brain", "enteric",
    "bowel", "ileum", "jejunum", "duodenum", "cecum", "caecum",
    "fecal", "faecal", "feces", "faeces", "probiotic", "prebiotic",
    "dysbiosis", "short chain fatty acid", "scfa", "leaky gut",
    "intestinal permeability", "epithelial barrier",
    "intestinal mucosa", "colonic mucosa", "duodenal mucosa", "rectal mucosa",
    "gut mucosa", "intestinal mucosal",
    "vagus nerve", "enteroendocrine",
]


GENE_EXPR_PATTERNS = [re.compile(p, re.IGNORECASE) for p in [
    r"\brna[\s\-]?seq(uencing)?\b",
    r"\bsc[\s\-]?rna[\s\-]?seq(uencing)?\b",
    r"\bsn[\s\-]?rna[\s\-]?seq(uencing)?\b",
    r"\bmrna[\s\-]?seq(uencing)?\b",
    r"\btotal[\s\-]rna[\s\-]?seq(uencing)?\b",
    r"\bbulk[\s\-]rna[\s\-]?seq(uencing)?\b",
    r"\blong[\s\-]read[\s\-]rna[\s\-]?seq\b",
    r"\brna\s+sequencing\b",
    r"\bsingle[\s\-]cell\s+rna[\s\-]?seq(uencing)?\b",
    r"\bsingle[\s\-]nucleus\s+rna[\s\-]?seq(uencing)?\b",
    r"\b(smart[\s\-]?seq\d*|drop[\s\-]?seq|mars[\s\-]?seq|cel[\s\-]?seq\d*|10[xX][\s\-]?(genomics|chromium)?|chromium[\s\-]?next\s*gem|parse[\s\-]?biosciences|split[\s\-]?pipe)\b",
    r"\bspatial[\s\-]?transcriptom(ics|e)\b",
    r"\b(visium|slide[\s\-]?seq|merfish|seqfish|stereo[\s\-]?seq|xenium)\b",
    r"\btranscriptom(ics|e|ic)\b",
    r"\bmicroarray\b",
    r"\bnanostring\b",
    r"\bgene\s+expression\b",
    r"\bdifferential(ly)?\s+(gene\s+)?expression\b",
    r"\bmrna\s+(expression|level|abundance)\b",
    r"\bdeseq\d*\b",
    r"\bedger\b",
    r"\blimma\b",
    r"\bvoom\b",
    r"\bseurat\b",
    r"\bscanpy\b",
    r"\bmonocle\b",
    r"\b(whole[\s\-]?)?transcriptome[\s\-]?sequencing\b",
    r"\bmeta[\s\-]?transcriptom(ics|ic|e)\b",
]]

MICROBIOME_PATTERNS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bmicrobiom(e|ics)\b",
    r"\bmicrobiota\b",
    r"\bgut[\s\-]?flora\b",
    r"\b16s[\s\-]?(r)?rna[\s\-]?(gene\s+)?sequencing\b",
    r"\b16s[\s\-]?(r)?rrna\b",
    r"\b16s[\s\-]?ribosomal\b",
    r"\bv[34][\s\-]?v[45]\s+(region|hypervariable)\b",
    r"\bits\d?\s+sequencing\b",
    r"\bamplicon[\s\-]?seq(uencing)?\b",
    r"\bshotgun[\s\-]?(metagenomic\s+)?sequencing\b",
    r"\bwhole[\s\-]?(metagen|genom)e[\s\-]?sequencing\b",
    r"\bwgs\b",
    r"\bmeta[\s\-]?genom(ics|ic|e)\b",
    r"\bmeta[\s\-]?transcriptom(ics|ic|e)\b",
    r"\b(operational\s+taxonomic\s+unit|otu)s?\b",
    r"\b(amplicon\s+sequence\s+variant|asv)s?\b",
    r"\b(alpha|beta|microbial|bacterial)\s+diversity\b",
    r"\b(shannon|simpson|chao1|faith|unifrac|bray[\s\-]?curtis)\b",
    r"\b(fecal|faecal|stool)\s+(microbiome|microbiota|transplant|transfer|16s)\b",
    r"\bfmt\b",
    r"\b(gut|intestinal|cecal|colonic)\s+(bacteria|microbi(al|ome|ota)|flora)\b",
    r"\bbacterial\s+composition\b",
    r"\bmicrobial\s+(composition|community|abundance|profile|diversity|metabolite)\b",
    r"\b(firmicutes|bacteroidetes|actinobacteria|proteobacteria|verrucomicrobia)\b",
    r"\b(lactobacillus|bifidobacterium|akkermansia|faecalibacterium|clostridium|prevotella|ruminococcus|bacteroides|blautia)\b",
    r"\b(short[\s\-]?chain\s+fatty\s+acid|scfa)s?\b",
    r"\b(butyrate|propionate|acetate|valerate)\b",
    r"\blipopolysaccharide\b",
    r"\blps[\s\-]?(induced|mediated|signaling|level)?\b",
    r"\bdysbiosis\b",
    r"\bleaky[\s\-]?gut\b",
    r"\b(intestinal|gut)\s+permeability\b",
    r"\b(probiotic|prebiotic|synbiotic)s?\b",
    r"\btryptophan[\s\-]?(metabolism|catabolism|pathway)\b",
    r"\bserotonin[\s\-]?(biosynthesis|synthesis|production)\b",
    r"\bkynurenine\s+(pathway|metabolism)\b",
]]

GENE_CONTEXT_PATTERNS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bexpress\w*\b",                           # expression, expressed, expressing, overexpression
    r"\bover[\s\-]?express\w*\b",                # overexpressed, overexpression
    r"\bdown[\s\-]?regulat\w*\b",                # downregulated, downregulation, down-regulate
    r"\bup[\s\-]?regulat\w*\b",                  # upregulated, upregulation, up-regulate
    r"\bdys[\s\-]?regulat\w*\b",                 # dysregulated, dysregulation
    r"\bregulat\w*\b",                           # regulated, regulation, regulator, regulatory
    r"\bgene\s+expression\b",
    r"\bdifferential(ly)?\s+(gene\s+)?expression\b",
    r"\bmrna\s+(expression|level|abundance)\b",
    # Protein / serum / plasma level measurements (e.g., ELISA, serum protein quantification)
    r"\blevel[s]?\s+(?:(?:of\s+)?\w[\w\s,\-/]{0,40}?\s+)?(?:were|was|are|is)\s+(?:significantly\s+)?(?:increased|decreased|elevated|reduced|higher|lower|altered)\b",
    r"\b(?:significantly\s+)?(?:increased|decreased|elevated|reduced|higher|lower)\s+(?:\w+\s+)?(?:serum|plasma|tissue|blood|protein|colonic|intestinal)\s+level[s]?\b",
    r"\bconcentration[s]?\s+(?:of\s+)?\w[\w\s,\-/]{0,30}?\s+(?:were|was)\s+(?:significantly\s+)?(?:increased|decreased|elevated|reduced|higher|lower)\b",
    r"\bELISA\b",
]]

_DRUG_KW = (
    r"drug[s]?|medication[s]?|antidepressant[s]?|antipsychotic[s]?|"
    r"anxiolytic[s]?|stimulant[s]?|vehicle|placebo|"
    # Antipsychotics
    r"haloperidol|clozapine|risperidone|olanzapine|quetiapine|"
    r"aripiprazole|ziprasidone|amisulpride|paliperidone|lurasidone|"
    r"chlorpromazine|fluphenazine|perphenazine|thioridazine|"
    # Antidepressants
    r"fluoxetine|sertraline|paroxetine|escitalopram|citalopram|"
    r"venlafaxine|duloxetine|mirtazapine|bupropion|trazodone|"
    r"imipramine|amitriptyline|nortriptyline|clomipramine|"
    # Mood stabilisers / anxiolytics / others
    r"lithium|licl|lithium\s+chloride|"
    r"valproate|valproic\s+acid|\bvpa\b|"
    r"lamotrigine|carbamazepine|\bcbz\b|"
    r"diazepam|lorazepam|midazolam|alprazolam|clonazepam|"
    r"ketamine|phencyclidine|pcp|"
    r"amphetamine|methylphenidate|ritalin|cocaine|morphine|heroin|"
    r"nicotine|cotinine|ethanol|alcohol|caffeine|"
    r"dexamethasone|corticosterone|glucocorticoid[s]?|hydrocortisone|"
    r"naloxone|naltrexone|"
    r"thapsigargin|rapamycin|cyclosporine|tacrolimus|"
    r"taltirelin|thymosin|cerebrolysin|"
    r"levomilnacipran|milnacipran|desvenlafaxine|reboxetine|"
    r"bcg|bacillus\s+calmette|adjuvant\s+therapy|"
    r"japanese\s+encephalitis\s+virus|\bjev\b|"
    r"bisphenol[\s\-]?a|\bbpa\b|"
    r"bisphenol[\s\-]?s|\bbps\b|"
    r"neuroleptic[s]?|"
    r"n[\s\-]?acetylcysteine|\bnac\b|"
    r"polyphenol[s]?|theogallin|resveratrol|curcumin|quercetin|"
    r"sigma[\s\-]?\d+\s+receptor\s+(agonist|antagonist)|"
    r"\bPRE084\b|\bDTG\b|\bBD1047\b|"
    r"herbal\s+(medicine|extract|formula|remedy|preparation)|"
    r"traditional\s+chinese\s+medicine|\bTCM\b|"
    r"oral\s+liquid|oral\s+decoction|"
    r"\bintragastric\s+administration\b|\bgavage\s+(?:administration|group)\b|"
    r"\btreatment\s+group\s+(?:received|was\s+given|were\s+given|was\s+administered)\b|"
    r"\bPills?\s+(?:were|was)\s+given\b|"
    r"plant\s+extract|phytochem\w+|"
    r"aloe[\s\-]vera|"
    r"ascorbic\s+acid|"
    r"propionic\s+acid|\bshort[\s\-]chain\s+fatty\s+acid[s]?\b|"
    r"wheat\s+(malt|germ)\s+extract|food[\s\-]derived\s+(source|extract)|"
    r"nutritional\s+strateg\w+|dietary\s+supplement\w*|"
    # Alzheimer's drugs
    r"memantine|donepezil|galantamine|rivastigmine|aducanumab|lecanemab|"
    # Parkinson's drugs
    r"levodopa|l[\s\-]?dopa|carbidopa|pramipexole|ropinirole|selegiline|rasagiline|"
    # Differentiation agents
    r"retinoic\s+acid|\bngf\s+treatment\b|"
    # Galactose model inducer / oxidative stress inducers
    r"d[\s\-]?galactose|hydrogen\s+peroxide|\bH2O2\b|"
    # Other common lab compounds used as interventions
    r"ceftriaxone|nanoparticl\w+|"
    r"vitamin\s+[ABCDE]\s+(supplement|supplementation|treatment|deficiency\s+rescue)"
)

# ── Model-inducer pool used by the abstract-level filter ──────────────────────
# If any of these disease-model inducers were administered to animals in the
# abstract, the whole abstract is excluded from gene-context scoring — the gene
# changes are a result of the experimental compound, not the disorder itself.
_MODEL_INDUCER_KW = (
    r"poly[\s\-]?i[\s\-:]?c|polyinosinic|polyriboinosinic[\s\-]polyribocytidilic|"
    r"mptp|1[\s\-]methyl[\s\-]4[\s\-]phenyl[\s\-]1,2,3,6[\s\-]tetrahydropyridine|"
    r"6[\s\-]ohda|6[\s\-]hydroxydopamine|"
    r"mk[\s\-]?801|dizocilpine|"
    r"rotenone|"
    r"streptozotocin|\bstz\b|"
    r"kainic\s+acid|kainate|"
    r"scopolamine|"
    r"reserpine|"
    r"pentylenetetrazol|\bptz\b|"
    r"cuprizone|"
    r"ibotenic\s+acid|"
    r"quinolinic\s+acid|"
    r"chlorpromazine|"
    r"pilocarpine|"
    r"lactacystin|"
    r"okadaic\s+acid|"
    r"aluminium\s+chloride|aluminum\s+chloride|\bacl3\b|"
    r"\blps\b|lipopolysaccharide|"
    r"bso|l[\s\-]buthionine|sulfoximine|"
    r"gbr[\s\-]?12909|"
    r"mk[\s\-]?801|dizocilpine|"
    r"trimethyltin|tmt|"
    r"colchicine|"
    r"paraquat|maneb"
)

# Matches "model-inducer ... administered/injected/treated" (either order, ≤120 chars apart)
_MI_ADMIN_RE = re.compile(
    r"\b(?:" + _MODEL_INDUCER_KW + r")\b.{0,120}"
    r"(?:treated|injected|administered|infused|exposed|received|given|induc\w+|challeng\w+)",
    re.IGNORECASE | re.DOTALL,
)
_MI_ADMIN_RE2 = re.compile(
    r"(?:treated|injected|administered|infused|exposed|received|given|challeng\w+).{0,120}"
    r"\b(?:" + _MODEL_INDUCER_KW + r")\b",
    re.IGNORECASE | re.DOTALL,
)
# Catches "X-treated", "X-induced", "X-challenged" hyphenated forms
_MI_HYPHEN_RE = re.compile(
    r"\b(?:" + _MODEL_INDUCER_KW + r")[\s\-](?:treated|induced|challenged)\b",
    re.IGNORECASE,
)
# Catches bare model inducer mentioned alongside any induction/stimulation language
_MI_BARE_RE = re.compile(
    r"\b(?:" + _MODEL_INDUCER_KW + r")\b.{0,200}"
    r"(?:stimulat\w+|induc\w+|trigger\w+|promot\w+|caus\w+)\s+\w[\w\s]{0,30}"
    r"(?:gene|mRNA|expression|transcript)",
    re.IGNORECASE | re.DOTALL,
)

def _abstract_uses_model_inducer(text):
    """Return True if the abstract describes animals administered a model-inducing compound."""
    return bool(
        _MI_ADMIN_RE.search(text)
        or _MI_ADMIN_RE2.search(text)
        or _MI_HYPHEN_RE.search(text)
        or _MI_BARE_RE.search(text)
    )

# Abstract-level check for therapeutic drug administration
_TD_ADMIN_VERBS = (
    r"treated|injected|administered|administration|infused|infusion|"
    r"exposed|exposure|received|given|pretreated|pretreatment|"
    r"co[\s\-]treated|prescribed|prescribing|underwent|undergone|applied"
)
_TD_ADMIN_RE = re.compile(
    r"\b(?:" + _DRUG_KW + r")\b.{0,150}(?:" + _TD_ADMIN_VERBS + r")",
    re.IGNORECASE | re.DOTALL,
)
_TD_ADMIN_RE2 = re.compile(
    r"(?:" + _TD_ADMIN_VERBS + r").{0,150}\b(?:" + _DRUG_KW + r")\b",
    re.IGNORECASE | re.DOTALL,
)
# Dosage units anywhere in the abstract = pharmacological study
_DOSAGE_RE = re.compile(
    r"\b\d+[\s.]?\d*[\s]*(mg/kg|mg\s*/\s*kg|μg|µg|μM|µM|nM|nmol|mg/ml|ug/kg|ng/kg)\b",
    re.IGNORECASE,
)

_TD_DRUG_ANYWHERE_RE  = re.compile(r"\b(?:" + _DRUG_KW + r")\b", re.IGNORECASE)
_TD_ADMIN_ANYWHERE_RE = re.compile(
    r"\b(?:injected|administered|administration|infused|infusion|"
    r"instilled|gavaged|implanted|pretreated|pretreatment|"
    r"challenged\s+with|challenge\s+with|"
    r"training|immunization|vaccination|inoculation)\b",
    re.IGNORECASE,
)
# Also catches "X-treated" / "X-induced" for any named drug (e.g. "LiCl-induced")
_TD_HYPHEN_RE = re.compile(
    r"\b(?:" + _DRUG_KW + r")[\s\-](?:treated|induced|mediated)\b",
    re.IGNORECASE,
)

def _abstract_uses_therapeutic_drug(text):
    """Return True if the abstract is a pharmacological/treatment study."""
    if _TD_ADMIN_RE.search(text) or _TD_ADMIN_RE2.search(text):
        return True
    if _DOSAGE_RE.search(text):
        return True
    if _TD_HYPHEN_RE.search(text):
        return True
    # Fallback: named drug mentioned anywhere + clear injection/admin verb anywhere
    if _TD_DRUG_ANYWHERE_RE.search(text) and _TD_ADMIN_ANYWHERE_RE.search(text):
        return True
    return False

# Non-pharmacological interventions — gene changes from these are also excluded
_NONDRUG_INTERVENTION_RE = re.compile(
    # ── Physical / neuromodulation procedures ─────────────────────────────────
    r"\b(electro)?acupuncture\b|"
    r"\belectroconvulsive\b|electric(al)?\s+(shock|stimulation)|"
    r"\bdeep\s+brain\s+stimulation\b|\bdbs\b|"
    r"\btranscranial\b|\btms\b|\btdcs\b|"
    r"\bphotobiomodulation\b|acousto[\s\-]optical|\bgamma\s+entrainment\b|"
    r"\b\d+[\s\-]hz\s+(stimulation|entrainment|gamma)\b|"
    r"\bvagal?\s+nerve\s+stimulation\b|"
    r"\blight\s+(?:therapy|treatment)\b|\bphototherapy\b|"
    r"\bultrasound\s+(stimulation|treatment|therapy)\b|"
    r"\btranscranial\s+ultrasound\b|"
    # ── Optogenetic / chemogenetic / viral gene manipulation ───────────────────
    r"\boptogeneti[ck]\w*\b|"
    r"\bchemogeneti[ck]\w*\b|\bdreadd\b|\bhm3d\b|\bhm4d\b|\b[ck]no\b|"
    r"\b(aav|lentivir\w+|adeno[\s\-]associated)\b|"
    r"\bviral\s+vector\b|\bgene\s+transfer\b|"
    r"\bsiRNAs?\b|\bshRNAs?\b|\bmiRNA\s+(mimic|inhibitor|injection)\b|"
    r"\bRNA\s+interference\b|\bRNAi\b|"
    r"\bCRISPR\b|\bcas9\b|"
    r"\bknockdown\s+(of|in|using|via|with|by)\b|"  # deliberate knockdown, not descriptive
    r"\boverexpression\s+vector\b|"
    # iPSC reprogramming/derivation methodology papers (not patient-derived disease studies)
    r"\b(?:reprogramming|reprogrammed)\s+(?:somatic\s+)?cells?\b|"
    r"\b(?:cord\s+blood|urine|saliva|hair\s+follicle)\s+(?:cells?\s+)?(?:for|as)\s+(?:a\s+)?(?:source|reprogramming)\b|"
    r"\biPSC\s+(?:derivation|generation|reprogramming|technology)\b|"
    r"\bgeneration\s+of\s+iPSC\b|"
    # In-vitro treatment: "treated the cells/neurons with X"
    r"\btreated\s+(?:the\s+)?(?:cells?|neurons?|cultures?|slices?|astrocytes?)\s+with\b|"
    r"\bstereotaxic\s+(injection|surgery|implant)\b|"
    r"\bintracranial\s+(injection|infusion|implant)\b|"
    r"\bretrograde\s+(tracing|virus|vector)\b|"
    r"\banterograde\s+(tracing|virus|vector)\b|"
    # ── Fear/learning paradigms (not disease models) ─────────────────────────
    r"\bfear\s+conditioning\b|\bconditioned\s+fear\b|"
    r"\bassociative\s+learning\b|\bpavlovian\b|"
    # ── Dietary / metabolic interventions ─────────────────────────────────────
    r"\bhigh[\s\-]fat\s+diet\b|\bketogenic\s+diet\b|\bcaloric\s+restriction\b|"
    r"\bdietary\s+(intervention|restriction|supplementation|deprivation)\b|"
    r"\bdietary\s+(pufa|omega[\s\-]?3|fatty\s+acid)\s+deprivation\b|"
    r"\bn[\s\-]?3\s+pufa\s+deprivation\b|"
    r"\bretraction\s+in\b|"  # retracted papers
    r"\bembryonic\s+stem\s+cell[s]?\b|"
    r"\bprobiotic[s]?\s*(supplementation|treatment|administration|intervention)?\b|"
    r"\bprebiotic[s]?\s*(supplementation|treatment|intervention)?\b|"
    r"\bsynbiotic[s]?\b|"
    r"\bintermittent\s+hypoxia\b|\bhypoxic\s+(pre)?conditioning\b|"
    r"\bhyperbaric\s+oxygen\b|"
    r"\bhyperoxia\b|\bneonatal\s+(hyperoxia|oxidative\s+stress|oxygen\s+exposure)\b|"
    r"\b\d+\s*%\s+oxygen\s+(exposure|treatment|for)\b|"
    r"\boxygen\s+exposure\b|"
    # ── Transplantation ───────────────────────────────────────────────────────
    r"\bfecal\s+(transplant|transfer|microbiota\s+transplant)\b|\bfmt\b|"
    r"\bbone\s+marrow\s+transplant\b|"
    # ── Surgery / implants ────────────────────────────────────────────────────
    r"\bcraniotomy\b|\blesion(ing)?\s+of\b|"
    r"\bvagotomy\b|\bcolectomy\b|"
    r"\btraumatic\s+brain\s+injury\b|\btbi\b|"
    r"\bcerebral\s+ischemia\b|\brain\s+ischemia\b|"
    r"\bstroke\s+model\b|\bmiddle\s+cerebral\s+artery\s+occlusion\b|\bmcao\b|"
    r"\bcommon\s+carotid\s+artery\s+occlusion\b|"
    # ── Exercise / physical training ─────────────────────────────────────────
    r"\btreadmill\b|\bvoluntary\s+(wheel\s+)?running\b|"
    r"\bexercise\s+(training|protocol|intervention|group)\b|"
    r"\benvironmental\s+enrichment\b|\benriched\s+environment\b|"
    # ── Cell / gene therapy ───────────────────────────────────────────────────
    r"\bstem\s+cell\s+(therapy|transplant|replacement|treatment|injection|infusion)\b|"
    r"\bcell\s+(replacement|transplant)\s+therapy\b|"
    r"\b(mesenchymal|pluripotent|neural|Muse)\s+(stem\s+)?cells?\s+(were|transplanted|injected|"
    r"administered|used\s+to|applied|is\s+a|are\s+a)\b|"
    r"\bMuse\s+cells?\b|"
    r"\bgene\s+therapy\b|"
    r"\bCAR[\s\-]T\b|"
    # ── Transgenic overexpression / knock-in models used as tools ─────────────
    r"\btransgenic\s+(livestock|animal|overexpression|construct)\b|"
    r"\bsite[\s\-]specific\s+integration\b|"
    r"\btransgene\s+(integration|expression|insert)\b|"
    # ── Irradiation ───────────────────────────────────────────────────────────
    r"\birradiation\b|\bradiation\s+(exposure|treatment)\b",
    re.IGNORECASE,
)

def _abstract_uses_nondrug_intervention(text):
    """Return True if the abstract studies gene changes caused by a non-drug intervention."""
    return bool(_NONDRUG_INTERVENTION_RE.search(text))

# ── Positive disease-model requirement ────────────────────────────────────────
# At least one of the following must be present for an abstract to pass.
# Covers: human patient tissue, postmortem, clinical samples, animal disease
# models (genetic or chemically-induced), and patient-derived iPSC/neurons.
_DISEASE_MODEL_RE = re.compile(
    # Human patients / clinical samples
    r"\bpatients?\s+with\b|\bdiagnosed\s+with\b|"
    r"\bsubjects?\s+with\b|\bindividuals?\s+with\b|"
    r"\bpostmortem\b|\bpost[\s\-]mortem\b|\bautopsy\b|\bbrain\s+bank\b|"
    r"\bbrain\s+tissue\s+from\b|\btissue\s+samples?\s+from\b|"
    r"\bpatient[\s\-](?:derived|specific)\b|"
    r"\bfrom\s+(?:patients?|affected\s+individuals?|diagnosed\s+(?:individuals?|subjects?))\b|"
    r"\bcase[\s\-]control\b|\bclinical\s+sample[s]?\b|"
    r"\bschizophrenia\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bbipolar\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bautism\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bADHD\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bdepression\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bdepressed\s+(?:patients?|subjects?|group|rats?|mice|mouse|animals?)\b|"
    r"\banxiety\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bPTSD\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bOCD\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bAlzheimer\w*\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    r"\bParkinson\w*\s+(?:patients?|subjects?|group|brain|cohort|sample|individual)\b|"
    # Abbreviation-first patterns (MDD, BPD, ASD, etc.)
    r"\bMDD\s+(?:patients?|subjects?|group|cohort|rats?|mice|model)\b|"
    r"\bBPD\s+(?:patients?|subjects?|group|cohort)\b|"
    r"\bASD\s+(?:patients?|subjects?|group|cohort|children|individuals?)\b|"
    # "patients with [disorder]" — covers any word order
    r"\bpatients?\s+with\s+(?:major\s+)?(?:depressive\s+disorder|depression|MDD|bipolar|schizophrenia|autism|ADHD|anxiety|PTSD|OCD|Alzheimer|Parkinson)\b|"
    r"\b(?:depressive|anxiety|autistic|schizophrenic)\s+(?:patients?|subjects?|individuals?|group|cohort)\b|"
    # Animal disease models — generic
    r"\b(?:mouse|rat|animal|primate|rodent|murine)\s+model\s+of\b|"
    r"\bmodel\s+(?:of|for)\s+\w[\w\s]{0,30}(?:disorder|disease|syndrome)\b|"
    r"\b(?:disorder|disease)[\s\-]model\b|"
    # Genetic animal models
    r"\b(?:knock[\s\-]?out|knock[\s\-]?in)\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\bKO\s+(?:mice|rats|mouse|rat)\b|"
    r"\btransgenic\s+(?:mice|rats|mouse|rat|animals?|model|line)\b|"
    r"\bmutant\s+(?:mice|rats|mouse|rat)\b|"
    r"\bhaploinsuffici\w+\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\b\w+[\s\-]deficient\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\b\w+[\s\-]knockout\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\b\w+\s+(?:KO|cKO|heterozygous|homozygous|null)\s+(?:mice|rats|mouse|rat)\b|"
    r"\b(?:conditional|inducible|tissue[\s\-]specific)\s+knockout\b|"
    r"\bgene[\s\-]targeted\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\bdisruption\s+of\s+the\s+\w+\s+gene\b|"
    r"\btargeted\s+disruption\s+of\b|"
    # Named genetic models
    r"\bAPP[\s\-/]?PS\d\b|\b5xFAD\b|\b3xTg[\s\-]?AD\b|"
    r"\bLRRK2\s+(?:G2019S|R1441|I2020T|Y1699C)\b|"
    r"\bShank\d[\s\-](?:KO|knockout|deficient|mutant)\b|"
    # Chemical disease model inducers (valid — animal shows disorder symptoms)
    r"\b6[\s\-]?(?:ohda|hydroxydopamine)\b|\bmptp\b|"
    r"\bstreptozotocin\b|\bstz\s+(?:mice|rats|model)\b|"
    r"\bpoly[\s\-]?i[\s\-:]?c\b|\bmk[\s\-]?801\b|"
    r"\bscopolamine\b|\bcuprizone\b|\bkainic\s+acid\b|"
    r"\bquinolinic\s+acid\b|\brotenone\b|\bparaquat\b|\breserpine\b|"
    r"\bokadaic\s+acid\b|\blactacystin\b|"
    # Patient-derived iPSC / neurons
    r"\biPSC[\s\-]?derived\b|\bpatient[\s\-]derived\s+(?:neurons?|cells?|iPSC)\b|"
    r"\bpatient[\s\-]specific\s+(?:neurons?|cells?|iPSC)\b|"
    r"\biPSC\s+from\s+(?:patients?|individuals?|subjects?)\b|"
    r"\bdisease[\s\-]specific\s+(?:neurons?|cells?|iPSC)\b|"
    # Stress-based animal models — valid for depression, anxiety, PTSD
    r"\bchronic\s+(?:mild|unpredictable|variable|restraint|social|defeat|immobilization)\s+stress\b|"
    r"\b(?:CMS|CUS|CUMS|CSDS)\s+(?:model|rats?|mice|mouse|group|protocol|paradigm|induced|exposed)\b|"
    r"\b(?:subjected|exposed)\s+to\s+(?:CMS|CUS|CUMS|CSDS|chronic\s+stress)\b|"
    r"\bchronic\s+stress\s+(?:model|rats?|mice|mouse|group|protocol|paradigm|induced)\b|"
    r"\bforced\s+swim\s+(?:test|model|stress)\b|\btail\s+suspension\s+(?:test|model)\b|"
    r"\brestraint\s+stress\s+(?:model|group|protocol|paradigm)\b|"
    r"\bsocial\s+defeat\s+(?:stress|model)\b|\bpredator\s+stress\s+(?:model|exposure)\b|"
    r"\bmaternal\s+separation\s+(?:model|stress|protocol)\b|"
    r"\bneonatal\s+(?:stress|separation|isolation)\s+(?:model|group|protocol)\b|"
    r"\bearly\s+life\s+(?:stress|adversity)\s+(?:model|group|protocol)\b|"
    r"\bimmobilization\s+stress\s+(?:model|group|protocol)\b",
    re.IGNORECASE,
)

def _has_disease_model(text):
    """Return True if the abstract contains an actual disease model (patient tissue or animal model)."""
    return bool(_DISEASE_MODEL_RE.search(text))

# Behavioral-primary exclusion — excludes abstracts whose explicit aim or primary
# result is behavioral/cognitive outcome, with gene expression only secondary.
_BEHAVIORAL_AIM_RE = re.compile(
    # Aim statement targeting behavior
    r"\baim(?:ed|s)?\s+(?:to\s+)?(?:was\s+to\s+)?(?:investigate|examine|assess|evaluate|study|"
    r"explore|determine|test|characterize|measure|quantify)\s+.{0,80}"
    r"(?:behavioral?|locomotor|cognitive|social\s+(?:behavior|interaction|memory)|"
    r"memory|learning|anxiety[\s\-]like|depression[\s\-]like|fear|attention)\s+"
    r"(?:deficits?|impairment|performance|phenotype|function|outcome|changes?|capacit\w+|abilit\w+)\b|"
    # "we investigated behavioral..." or "we assessed cognitive..."
    r"\bwe\s+(?:investigated|examined|assessed|evaluated|studied|explored|measured|tested|"
    r"characterized)\s+.{0,60}"
    r"(?:behavioral?|locomotor|cognitive|social)\s+"
    r"(?:deficits?|impairment|performance|phenotype|function|outcome|changes?)\b|"
    # Results/conclusion centred on behavioral outcomes
    r"\b(?:behavioral?|locomotor|cognitive|social)\s+"
    r"(?:deficits?|impairment|performance|changes?|phenotype)\s+"
    r"(?:were|was)\s+(?:the\s+)?(?:primary|main|principal|key|major)\s+"
    r"(?:endpoint|outcome|finding|result|readout)\b|"
    # "behavioral deficits were observed / demonstrated / found" as a conclusion statement
    # combined with gene expression mentioned only in passing (≤1 expression-related word nearby)
    r"\bbehavioral?\s+(?:deficits?|impairment|changes?|outcomes?)\s+were\s+"
    r"(?:observed|found|detected|demonstrated|shown|noted|documented)\b",
    re.IGNORECASE | re.DOTALL,
)

# Require at least one explicit gene expression measurement signal
_EXPRESSION_MEASUREMENT_RE = re.compile(
    # Measurement methods
    r"\bmRNA\b|"
    r"\bRT[\s\-]?PCR\b|\bqPCR\b|q[\s\-]?RT[\s\-]?PCR\b|"
    r"\bmicroarray\b|\bRNA[\s\-]?seq\b|\bscRNA[\s\-]?seq\b|"
    r"\btranscriptom\w+\b|"
    r"\bin\s+situ\s+hybridization\b|\bISH\b|"
    r"\bWestern\s+blot\b|\bimmunoblot\b|"
    r"\bimmunohistochem\w+\b|\bimmunofluo\w+\b|"
    r"\bexpression\s+profil\w+\b|"
    # Gene expression change descriptors — verb forms
    r"\b(?:up|down)[\s\-]?regulat\w+\b|"
    r"\bdifferentially\s+expressed\b|\bdifferential\s+(?:gene\s+)?expression\b|"
    r"\boverexpressed\b|\bunderexpressed\b|"
    r"\bexpression\s+(?:was|were|is|are)\s+\w*\s*"
    r"(?:increased|decreased|reduced|elevated|altered|changed|upregulated|"
    r"downregulated|measured|quantified|higher|lower|suppressed|enhanced|"
    r"significantly\s+\w+)\b|"
    # "increased/decreased/reduced/elevated/altered expression" (adjective before noun)
    r"\b(?:increased|decreased|reduced|elevated|altered|changed|higher|lower|"
    r"enhanced|suppressed|significant)\s+(?:\w+\s+){0,3}expression\b|"
    # "expression levels increased/decreased/..."
    r"\bexpression\s+levels?\s+(?:were|was|are|is)?\s*"
    r"(?:increased|decreased|reduced|elevated|altered|changed|higher|lower|"
    r"significantly\s+\w+)\b|"
    # "mRNA levels were..."
    r"\bmRNA\s+levels?\s+(?:were|was|are|is)?\s*"
    r"(?:increased|decreased|reduced|elevated|altered|changed|higher|lower)\b|"
    # gene expression analysis / study / data
    r"\bgene\s+expression\s+(?:analysis|profil\w+|data|levels?|changes?|pattern|studi\w+)\b|"
    # Protein quantification methods
    r"\bELISA\b|\bimmunoas{2}ay\b|\bimmuno[\s\-]?precipitation\b|"
    # Protein/serum/plasma levels measured — flexible for comma-separated lists
    # e.g. "plasma levels of zonulin, I-FABP, LPS, and claudin-5 were significantly higher"
    r"\b(?:protein|serum|plasma|tissue|colonic|intestinal|blood|fecal|urine|CSF)\s+level[s]?\s+"
    r"(?:of\s+)?.{0,80}?\s+(?:were|was|are|is)\s+(?:significantly\s+)?(?:increased|decreased|reduced|elevated|higher|lower|altered)\b|"
    r"\blevel[s]?\s+(?:of\s+)?.{0,60}?\s+(?:were|was)\s+(?:significantly\s+)?(?:increased|decreased|reduced|elevated|higher|lower)\b|"
    # "X levels were significantly Y" without the preceding qualifier
    r"\b\w[\w\-/]+\s+level[s]?\s+(?:were|was|are|is)\s+(?:significantly\s+)?(?:increased|decreased|reduced|elevated|higher|lower|altered)\b|"
    # "concentrations of X were"
    r"\bconcentration[s]?\s+(?:of\s+)?.{0,60}?\s+(?:were|was)\s+(?:significantly\s+)?(?:increased|decreased|elevated|reduced|higher|lower)\b|"
    # Reversed: "significantly lower/higher serum/plasma levels of X"
    r"\b(?:significantly\s+)?(?:increased|decreased|reduced|elevated|higher|lower|altered)\s+"
    r"(?:\w+\s+)?(?:serum|plasma|tissue|blood|protein|CSF|fecal|urinary)\s+level[s]?\b|"
    # Gut-specific measurement methods
    r"\bintestinal\s+permeability\b|\bgut\s+permeability\b|\bleaky[\s\-]?gut\b|"
    r"\bFITC[\s\-]?dextran\b|\blactulose[\s\-]?mannitol\b|"
    r"\btransepithelial\s+(?:electrical\s+)?resistance\b|\bTEER\b|"
    r"\bzonulin\s+(?:levels?|concentration|serum|plasma|fecal|stool)\b|"
    r"\bcalprotectin\b|\bI[\s\-]?FABP\b|\bfatty\s+acid[\s\-]?binding\s+protein\b|"
    r"\b(?:fecal|stool|intestinal|colonic|mucosal)\s+(?:levels?\s+of|concentration\s+of|abundance\s+of)\b|"
    r"\bmucosa[l]?\s+(?:integrity|damage|inflammation|biopsy|sample)\b|"
    r"\bcolonic\s+(?:biopsy|tissue|mucosa|epithelium|sample)\b|"
    r"\benteric\s+(?:nervous\s+system|neurons?|glial)\b|"
    r"\bgut[\s\-]?brain\s+(?:axis|interaction|communication|signaling)\b",
    re.IGNORECASE,
)

def _abstract_is_behavioral_primary(text):
    """Return True if abstract's primary aim/result is behavioral with gene expression secondary."""
    if not _BEHAVIORAL_AIM_RE.search(text):
        return False
    return not _EXPRESSION_MEASUREMENT_RE.search(text)

# Genetic association / polymorphism studies — primary outcome is disease risk or
# genotype–phenotype association, not gene expression levels.
_GENETIC_ASSOC_RE = re.compile(
    # Polymorphism/SNP/genotype as subject of association
    r"\b(?:polymorphism[s]?|SNP[s]?|genotype[s]?|haplotype[s]?|allele[s]?|variant[s]?|"
    r"copy\s+number\s+variation[s]?|\bCNV[s]?\b)\s+.{0,60}"
    r"(?:(?:is|are|was|were)\s+(?:not\s+)?(?:significantly\s+)?associated\s+with|"
    r"(?:significantly\s+)?predict[s]?\s+risk|confer[s]?\s+(?:increased\s+)?risk)\b|"
    # GWAS / linkage / association study language
    r"\bgenome[\s\-]wide\s+association\b|\bGWAS\b|"
    r"\blinkage\s+(?:disequilibrium|analysis|study)\b|\bLD\s+block\b|"
    r"\bassociation\s+(?:study|analysis|between\s+\w+\s+(?:polymorphism|SNP|genotype))\b|"
    r"\bno\s+(?:significant\s+)?association\s+(?:was\s+found\s+)?between\s+\w[\w\s]{0,30}"
    r"(?:polymorphism|SNP|genotype|variant|allele)\b|"
    # Risk allele / susceptibility gene framing (pure genetics)
    r"\brisk\s+(?:allele[s]?|variant[s]?|genotype[s]?|haplotype[s]?)\b|"
    r"\bsusceptibility\s+(?:allele[s]?|variant[s]?|locus|loci)\b|"
    r"\b(?:disease|disorder)\s+risk\s+(?:allele|variant|SNP|genotype)\b",
    re.IGNORECASE | re.DOTALL,
)

# Absolute exclusion: polymorphism/SNP as the FINDING (regardless of background mentions)
_POLYMORPHISM_FINDING_RE = re.compile(
    r"\b(?:polymorphism[s]?|SNP[s]?|genotype[s]?|haplotype[s]?|allele[s]?)\s+.{0,100}"
    r"(?:(?:is|are|was|were)\s+(?:not\s+)?(?:significantly\s+)?associated\s+with|"
    r"confer[s]?\s+(?:increased\s+)?risk|predict[s]?\s+(?:the\s+)?risk)\b|"
    r"\bno\s+(?:significant\s+)?association\s+(?:was\s+)?(?:found|detected|observed)?\s+between\b|"
    r"\bGWAS\b|\bgenome[\s\-]wide\s+association\b",
    re.IGNORECASE | re.DOTALL,
)

def _get_results_section(text):
    """Return the RESULTS/CONCLUSIONS portion of an abstract, or the latter 55% if unstructured."""
    m = re.search(
        r'\b(?:RESULTS?|FINDING[S]?|CONCLUSION[S]?|OUTCOME[S]?|WE\s+FOUND|WE\s+SHOW|WE\s+REPORT|'
        r'WE\s+DEMONSTRATE|WE\s+OBSERVED|WE\s+IDENTIFIED|IN\s+THIS\s+STUDY,?\s+WE)'
        r'[\s:]+',
        text, re.IGNORECASE,
    )
    if m:
        return text[m.start():]
    return text[int(len(text) * 0.45):]

def _is_genetic_association_study(text):
    """Return True if abstract's FINDINGS are about polymorphism/SNP association, not expression."""
    results = _get_results_section(text)
    # Polymorphism as a finding = absolute exclusion
    if _POLYMORPHISM_FINDING_RE.search(results):
        return True
    # Broader genetic association study without expression measurement anywhere
    if _GENETIC_ASSOC_RE.search(text) and not _EXPRESSION_MEASUREMENT_RE.search(text):
        return True
    # Polymorphism is the PRIMARY finding even if expression is mentioned in passing
    # e.g. "DRD2 polymorphism was consistently associated with gait function"
    _POLY_PRIMARY_RE = re.compile(
        r"\b(?:polymorphism[s]?|SNP[s]?|genotype[s]?|haplotype[s]?|allele[s]?)\s+.{0,150}"
        r"(?:(?:was|were|is|are)\s+(?:consistently|significantly|strongly|independently|positively|negatively|inversely)?\s*associated\s+with|"
        r"predict[s]?\s+|confer[s]?\s+)",
        re.IGNORECASE | re.DOTALL,
    )
    if _POLY_PRIMARY_RE.search(results):
        return True
    return False

# ── Computational / bioinformatics studies — no original measurements ─────────
# These use public databases or PPI network analysis without measuring themselves.
_COMPUTATIONAL_STUDY_RE = re.compile(
    r"\bppi\s+(?:network|sub[\s\-]?network)\b|"
    r"\bprotein[\s\-]protein\s+interaction\s+network\b|"
    r"\bco[\s\-]?expression\s+(?:network|analysis)\s+(?:based\s+on\s+)?public\b|"
    r"\bpublicly\s+available\s+(?:microarray|expression|gene\s+expression)\s+(?:data|dataset)\b|"
    r"\bpublic\s+(?:microarray|gene\s+expression|expression)\s+(?:database|data|dataset)\b|"
    r"\b(?:GEO\b|ArrayExpress|TCGA|GTEx)\s+(?:database|dataset|data|repository)\b|"
    r"\bin\s+silico\s+(?:analysis|study|prediction|screening)\b|"
    r"\bbioinformatics\s+(?:analysis|study|approach|pipeline|workflow)\b.{0,100}"
    r"(?:candidate\s+genes?|network|pathway)\b|"
    r"\bnearest\s+neighbour[s]?\s+of\s+(?:the\s+)?\d+\s+(?:reported\s+)?candidate\s+genes?\b",
    re.IGNORECASE | re.DOTALL,
)

def _is_computational_study(text):
    """Return True if the abstract is a computational/bioinformatics study using public data."""
    return bool(_COMPUTATIONAL_STUDY_RE.search(text))

# ── Fluid-only measurement — serum/plasma/CSF/blood/urine is not brain/gut tissue ──
# If the abstract ONLY measures the gene in body fluids (not tissue), we exclude it.
_FLUID_MEASUREMENT_RE = re.compile(
    r"\b(?:serum|plasma|blood|CSF|cerebrospinal\s+fluid|urine|urinary|saliva|salivary)\s+"
    r"(?:levels?\s+of|concentrations?\s+of|samples?\s+of|measurement\s+of|"
    r"(?:were|was)\s+(?:measured|quantified|assessed|analyzed))\b",
    re.IGNORECASE,
)
_TISSUE_MEASUREMENT_RE = re.compile(
    r"\b(?:brain|cortex|hippocampus|hippocampal|striatum|striatal|prefrontal|amygdala|"
    r"cerebellum|substantia\s+nigra|basal\s+ganglia|spinal\s+cord|"
    r"intestine|intestinal|colon|colonic|gut|ileum|jejunum|duodenum|"
    r"postmortem|post[\s\-]mortem|biopsy|tissue\s+sample|"
    r"neurons?|astrocytes?|microglia)\s+"
    r"(?:tissue|sample|biopsy|extract|lysate|homogenate|expression|levels?|"
    r"mRNA|protein|immunostaining|immunohistochem\w+)\b|"
    r"\b(?:tissue|biopsy|postmortem|post[\s\-]mortem)\s+(?:from|of|in)\b",
    re.IGNORECASE,
)

def _is_fluid_only_measurement(text):
    """Return True if gene expression is measured ONLY in body fluids, not brain/gut tissue."""
    if not _FLUID_MEASUREMENT_RE.search(text):
        return False
    # If there is also a tissue measurement, allow it through
    if _TISSUE_MEASUREMENT_RE.search(text):
        return False
    return True

# Blood-derived cell lines — not brain/gut tissue
_BLOOD_CELL_LINE_RE = re.compile(
    r"\blymphoblastoid\s+cell\s+(?:lines?|cultures?)\b|"
    r"\bLCLs?\b|"
    r"\bperipheral\s+blood\s+(?:mononuclear\s+cells?|lymphocytes?|leukocytes?|monocytes?|"
    r"granulocytes?|neutrophils?)\b|"
    r"\bPBMCs?\b|"
    r"\bwhole\s+blood\s+(?:RNA|samples?|cells?|transcriptom\w+|gene\s+expression)\b|"
    r"\bblood\s+(?:monocytes?|lymphocytes?|leukocytes?|neutrophils?|granulocytes?|"
    r"platelets?|mononuclear\s+cells?)\b",
    re.IGNORECASE,
)

def _is_blood_cell_line_study(text):
    """Return True if the abstract primarily uses blood-derived cell lines (not brain/gut tissue)."""
    if not _BLOOD_CELL_LINE_RE.search(text):
        return False
    if _TISSUE_MEASUREMENT_RE.search(text):
        return False
    return True

# ── Long-term depression = LTD (synaptic plasticity), not Major Depression ────
_LTD_SYNAPTIC_RE = re.compile(
    r"\blong[\s\-]term\s+depression\b.{0,300}\b(?:LTP|long[\s\-]term\s+potentiation|synaptic\s+plasticity|"
    r"hippocampal\s+(?:LTP|plasticity)|dendritic|NMDA|AMPA|glutamatergic)\b|"
    r"\b(?:LTP|long[\s\-]term\s+potentiation|synaptic\s+plasticity|hippocampal\s+plasticity|"
    r"dendritic\s+spines?|NMDA|AMPA)\b.{0,300}\blong[\s\-]term\s+depression\b",
    re.IGNORECASE | re.DOTALL,
)

def _abstract_is_ltd_not_mdd(text, disorder):
    """For Major Depression abstracts: return True if 'long-term depression' means LTD, not MDD."""
    if disorder != "Major Depression":
        return False
    # Only a problem if long-term depression appears and no MDD-specific language exists
    if not re.search(r"\blong[\s\-]term\s+depression\b", text, re.IGNORECASE):
        return False
    if re.search(r"\b(?:major\s+depressive\s+disorder|MDD|depressed\s+patients?|"
                 r"depression\s+(?:patients?|group|subjects?|scores?|symptoms?|diagnosis))\b",
                 text, re.IGNORECASE):
        return False  # genuine MDD paper with LTD also mentioned
    return bool(_LTD_SYNAPTIC_RE.search(text))

# Epigenetics / miRNA exclusion — this study is about gene expression only
_EPIGENETICS_RE = re.compile(
    r"\bepigeneti[ck]\w*\b|\bepigenome\b|\bepigenomic\w*\b|"
    # methylation and variants: methylation, demethylation, methylated, unmethylated,
    # hypermethylation, hypomethylation, remethylation etc.
    r"\b(de|un|re|hyper|hypo)?methylat\w+\b|"
    r"\bDNA\s+methylation\b|\bCpG\b|"
    # acetylation variants: acetylation, deacetylation, acetylated, deacetylated
    r"\b(de)?acetylat\w+\b|"
    # ubiquitination variants
    r"\bubiquitinat\w+\b|\bubiquitylat\w+\b|"
    # histone marks / chromatin
    r"\bhistone\s+(modification|mark|acetyl|methyl|ubiquitin|variant)\w*\b|"
    r"\bchromatin\s+(remodel\w+|modifi\w+|accessib\w+|structure)\b|"
    r"\bchromatin\s+immunoprecipitation\b|\bChIP\b|"
    # miRNA / small RNA / non-coding RNA
    r"\bmiRNA[s\-]?\b|\bmicro[\s\-]?RNA[s]?\b|\bmiR[\s\-]\d+\b|"
    r"\bsRNA[s]?\b|\bsmall[\s\-]RNA[s]?\b|"
    r"\bnon[\s\-]coding[\s\-]RNA[s]?\b|\bncRNA[s]?\b|"
    r"\blncRNA[s]?\b|\blong[\s\-]non[\s\-]coding[\s\-]RNA[s]?\b|"
    r"\bcircRNA[s]?\b|\bpiRNA[s]?\b|"
    r"\bSUMOylation\b|\bsumoylat\w+\b",
    re.IGNORECASE,
)

def _abstract_uses_epigenetics(text):
    """Return True if the abstract is about epigenetics or miRNA rather than gene expression."""
    return bool(_EPIGENETICS_RE.search(text))

# Catches treatment-attribution language regardless of what the compound is:
#   "downregulated by zymosan treatment"
#   "MLi-2-mediated gene expression"
#   "effects of [X] on gene/mRNA expression"
#   "X treatment altered/changed/affected expression"
# Named drug appearing next to "treatment" (noun) as a strong signal
_DRUG_TREATMENT_NOUN_RE = re.compile(
    r"\b(?:" + _DRUG_KW + r")\s+treatment\b|"
    r"\btreatment\s+with\s+(?:" + _DRUG_KW + r")\b|"
    r"\bduring\s+(?:" + _DRUG_KW + r")\s+treatment\b",
    re.IGNORECASE,
)
# "Li treatment" / "Li-treated" — lithium abbreviation
_LI_TREATMENT_RE = re.compile(
    r"\bLi[\s\-](treatment|treated|pretreatment|pretreated)\b",
    re.IGNORECASE,
)

_THERAPEUTIC_EFFECTS_RE = re.compile(
    r"\btherapeutic\s+effects?\s+of\s+(?:" + _DRUG_KW + r")\b|"
    r"\bsmall\s+molecule\s+\w[\w\-]+\s+(to\s+)?(activate|inhibit|target|treat|rescue)\b|"
    r"\b\w[\w\-]+,?\s+a\s+(small\s+molecule|TrkB\s+agonist|kinase\s+inhibitor|receptor\s+agonist"
    r"|receptor\s+antagonist|potent\s+and\s+specific)\b",
    re.IGNORECASE,
)

_TREATMENT_ATTRIBUTION_RE = re.compile(
    r"\b(up|down)[\s\-]?regulated\s+by\s+\w[\w\s,\-]+treatment\b|"
    r"\b\w[\w\-]+[\s\-]mediated\s+(gene\s+|mRNA\s+|differential\s+)?expression\b|"
    r"\b\w[\w\-]+[\s\-]induced\s+(changes?\s+in\s+)?(gene\s+|mRNA\s+|differential\s+)?expression\b|"
    r"\beffects?\s+of\s+\w[\w\s,\-]+\s+on\s+"
    r"(gene|mRNA|differential|protein|cellular|transcriptomic|genomic|"
    r"epigenetic|signaling|molecular|dopamine|serotonin)\b|"
    r"\b\w[\w\-]+\s+treatment\s+(significantly\s+)?(altered?|changed?|affected?|modulated?|"
    r"increased?|decreased?|upregulated?|downregulated?|regulated?|reversed?|restored?|"
    r"attenuated?|enhanced?|reduced?|elevated?|normalized?|rescued?|prevented?|improved?)\b|"
    # "treatment with X [change verb]" — catches unknown drugs too
    r"\btreatment\s+with\s+[\w\s,\-\+\(\)]{2,60}"
    r"(significantly\s+)?(increased?|decreased?|reversed?|improved?|restored?|attenuated?|"
    r"enhanced?|reduced?|elevated?|normalized?|rescued?|prevented?|altered?|modulated?)\b|"
    r"\b(treated|treatment)\s+with\s+\w[\w\s,\-]+(increased?|decreased?|altered?|"
    r"upregulated?|downregulated?|modulated?)\s+(gene\s+|mRNA\s+|protein\s+)?expression\b|"
    r"\bchallenged\s+with\s+\w[\w\s\-,]+(?:inhibitor|agonist|antagonist|compound|drug|agent)\b|"
    r"\btreatment\s+with\s+(the\s+)?recombinant\b|"
    r"\btreated\s+with\s+(the\s+)?recombinant\b|"
    r"\brecombinant\s+\w[\w\-]+\s+(protein|peptide|factor|cytokine)\b|"
    # drug directly causing expression change without 'treatment' word
    r"\b(?:" + _DRUG_KW + r")\s+(significantly\s+)?"
    r"(increased?|decreased?|upregulated?|downregulated?|altered?|modulated?|"
    r"induced?|suppressed?|enhanced?|reduced?)\s+(the\s+)?"
    r"(gene\s+|mRNA\s+|protein\s+)?expression\b",
    re.IGNORECASE | re.DOTALL,
)

def _abstract_attributes_expression_to_treatment(text):
    """Return True if gene expression changes are explicitly attributed to any compound/treatment."""
    return bool(
        _TREATMENT_ATTRIBUTION_RE.search(text)
        or _THERAPEUTIC_EFFECTS_RE.search(text)
        or _DRUG_TREATMENT_NOUN_RE.search(text)
        or _LI_TREATMENT_RE.search(text)
    )

# Sentences matching this are skipped — they describe gene changes CAUSED BY
# treatment/drugs, not by the disorder itself.
#
# Changes vs. original:
#   • "induced by" and "in response to" are now handled by helper functions
#     (_induced_by_drug / _in_response_to_drug) that require a drug keyword in
#     the same sentence, so "stress-induced" or "disease-induced" no longer match.
#   • Added explicit drug/medication keywords as causal triggers.
#   • Added standalone named drug patterns (haloperidol, clozapine, etc.) so
#     sentences that describe what a specific drug did are always excluded.
TREATMENT_CAUSAL_RE = re.compile(
    # Direct administration language
    r"\btreated\s+with\b"
    r"|\btreatment\s+with\b"
    r"|\badministration\s+of\b"
    # Reversal / rescue attributed to a treatment
    r"|\b(reversed|restored|normalized|attenuated|rescued|ameliorated)\s+by\b"
    # Drug dosage units — dead giveaway of pharmacological study
    r"|\b\d+[\s]*(mg/kg|μg|μM|nM|mg/ml|nmol|µg|µM)\b"
    # Injection / infusion protocols
    r"|\bfollowing\s+\w+\s*(treatment|administration|injection|infusion)\b"
    r"|\bafter\s+\w+\s*(treatment|administration|injection)\b"
    # Drug keyword DIRECTLY preceding a change verb (e.g. "haloperidol induced", "drugs increased")
    r"|(?:" + _DRUG_KW + r")\s*[\-\s]?\s*(induced|altered|changed|increased|decreased|elevated|reduced)\b"
    # Change verb attributed BY to a drug (e.g. "decreased by haloperidol", "altered by drugs")
    r"|\b(induced|altered|changed|increased|decreased|elevated|reduced)\s+by\s+(?:" + _DRUG_KW + r")\b",
    re.IGNORECASE
)

def _treatment_sentence(sent):
    """Return True if the sentence describes a drug/treatment causing the gene change.

    Handles the two cases the compiled regex intentionally excludes:
      • 'induced by'   → only drug-induced counts; stress/disease-induced does not.
      • 'in response to' → only when a drug keyword is in the same sentence.
    """
    if TREATMENT_CAUSAL_RE.search(sent):
        return True
    # "induced by <drug>" — but NOT "induced by stress" or "induced by the disease"
    if re.search(r"\binduced\s+by\b", sent, re.IGNORECASE) and \
            re.search(r"\b(?:" + _DRUG_KW + r")\b", sent, re.IGNORECASE):
        return True
    # "in response to <drug>" — but NOT "in response to chronic stress"
    if re.search(r"\bin\s+response\s+to\b", sent, re.IGNORECASE) and \
            re.search(r"\b(?:" + _DRUG_KW + r")\b", sent, re.IGNORECASE):
        return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# (2) GENE LISTS
#
# Gene detection: instead of a fixed candidate list, we detect any gene-like
# symbol that appears in a biological context sentence in a passing abstract.
# Tissue classification is purely context-based (brain/gut keywords near the
# gene mention), so no predefined lists are needed.
# ─────────────────────────────────────────────────────────────────────────────

# Regex to find potential gene symbols: uppercase, 2–9 chars, may contain digits
_GENE_SYMBOL_RE = re.compile(r'\b([A-Z][A-Z0-9]{1,8})\b')

# Non-gene abbreviations that match the pattern but are not gene symbols
_NON_GENE_SYMBOLS = frozenset({
    # Nucleic acids / molecular biology techniques
    "DNA","RNA","MRNA","RRNA","TRNA","SNRNA","LNCRNA","MIRNA","SIRNA","SHRNA",
    "NCRNA","CDNA","GDNA","MTDNA","PCR","RTPCR","QPCR","NGS","WGS","WES",
    "CHIP","ATAC","RNAI","CRISPR","CAS","CAS9","CHIP","CAGE","CLIP",
    # Techniques / assays
    "ELISA","IHC","ICC","IF","WB","IP","CO","FACS","FISH","ISH","EMSA",
    "NMR","MRI","FMRI","PET","CT","EEG","EMG","ECG","EKG","MEG","SPECT",
    "LC","LCMS","HPLC","GC","SDS","PAGE","SDSPAGE","MALDI","MSMS",
    "RT","IHC","TUNEL","BCA","MTT","LDH","RIA","CLIA",
    # Clinical / psychiatric disorder abbreviations
    "MDD","BPD","ASD","ADHD","OCD","GAD","SAD","PTSD","TBI","ALS","HD",
    "SCZ","SZP","BD","IBD","IBS","CD","UC","MS","RA","SLE","CFS",
    # Epidemiology / statistics
    "BMI","BP","HR","RR","OR","CI","SD","SE","SEM","IQR","AUC","ROC",
    "ANOVA","ANCOVA","MANOVA","GLM","LME","GEE","PCA","ICA","SVM","ANN",
    "RCT","ITT","FDR","FPR","SNP","CNV","GWAS","QTL","LOD","HWE",
    # Routes of administration / dosing
    "IV","IM","SC","PO","IP","SQ","IT","IO",
    # Units
    "MG","UG","NG","PG","ML","UL","NL","MM","UM","NM","PM","KG","LB",
    "MHZ","KHZ","KDA","GHZ","MS","US","NS","PS","MPH","BPM",
    # Biochemistry / metabolites
    "ATP","ADP","AMP","GTP","GDP","GMP","CAMP","CGMP","NAD","NADH","NADPH",
    "FAD","FADH","FMN","COA","ACOA","ROS","NO","CO2","H2O","O2","H2O2",
    "PBS","BSA","DTT","EDTA","HEPES","TRIS","DMSO","DMEM","RPMI","HBSS",
    "LPS","LTA","CFA","IFA","AGO","PHA","ConA",
    # Organizations / countries
    "FDA","WHO","NIH","CDC","EU","USA","UK","US","UN","NATO",
    # Neuroanatomy regions (abbreviated, not gene symbols)
    "CNS","PNS","BBB","CSF","ENS","ANS","SNS","PFC","ACC","OFC","DLPFC",
    "VTA","SN","SNc","SNr","NAC","STR","HIP","AMY","CTX","CB","BS",
    "DRN","LC","MRN","PVN","SON","BLA","CeA","BMA","HPA","SAM","HPG",
    "CA1","CA2","CA3","CA4","DG","EC","SC","IC","TC",
    # GI anatomy (not gene names)
    "GI","GIT","IBD","SIBO","GER","GERD","IEL","LP","PP","MLN","GLP",
    # Hormones / classes (not gene symbols)
    "HRT","OCP","SSRI","SNRI","MAOI","TCA","NRI","NDRI","SGA","FGA",
    # Abstract section headers
    "BACKGROUND","OBJECTIVE","OBJECTIVES","PURPOSE","INTRODUCTION","AIM","AIMS",
    "CONTEXT","SUMMARY","METHODS","METHOD","RESULTS","RESULT","CONCLUSIONS",
    "CONCLUSION","ABSTRACT","FINDINGS","DISCUSSION","SIGNIFICANCE","KEYWORDS",
    # Common abbreviations in running text
    "THE","AND","FOR","NOT","ARE","WAS","WERE","HAS","HAD","WITH","FROM",
    "THIS","THAT","THEY","THEM","BEEN","HAVE","INTO","THAN","THEN","WHEN",
    "BUT","CAN","ITS","OUR","ALL","MAY","NEW","TWO","ONE","ANY","HOW",
    "VS","ET","AL","IE","EG","IN","OF","ON","AT","BY","TO","AN","AS","IS",
    # Protein/gene classes (not individual gene symbols)
    "GPCR","RTK","NTR","VGCC","TRP","HCN","ERK","JNK","RAS","RAC","CDC",
    "CDK","CYC","PKA","PKC","PKG","PKB","PDE","PLC","PLA","DAG","IP","PKC",
})

# Kept for common-name → gene symbol mapping (e.g. "occludin" → OCLN)
# These are still used by the alias detection path.
BRAIN_GENES = [
    # ── Neuronal identity / pan-neuronal ──────────────────────────────────────
    "RBFOX3","MAP2","TUBB3","SYP","SYT1","DLG4","GFAP","MBP","AIF1","ENO2",
    "NEFH","NEFL","SNAP25","TH","GAD1","SLC17A7","PVALB","SST","OLIG2","S100B",
    "DCX","STMN2","ASCL1","PROX1","NRGN","CAMK2A","SLC17A6","GAD2","SLC32A1",
    "FOXP2","ALDH1L1","SLC1A2","AQP4","MOG","PLP1","CNP","TMEM119","P2RY12","CX3CR1",
    "MAG","MOBP","MYRF","OLIG1","GJA1","CD68","TREM2","ITGAM",
    # ── Neurotrophins & growth factors ────────────────────────────────────────
    "BDNF","NGF","NTF3","NTF4","NTRK1","NTRK2","NTRK3","GDNF","CNTF",
    "FGF2","FGF1","EGF","VEGFA","VEGFB","IGF1","IGF2","IGF1R","IGFBP3",
    "PPARGC1A","TFAM","SIRT1","SIRT3",
    # ── Neuropeptides & receptors ─────────────────────────────────────────────
    "NPY","VIP","GAL","ADCYAP1","ADCYAP1R1","CALCA","CALCB",
    "TAC1","PENK","PDYN","POMC","AGRP","CRH","CRHR1","CRHR2","CORT",
    "AVP","AVPR1A","AVPR1B","AVPR2","OXT","OXTR",
    "CCK","NTS","GRP","GRPR","GALR1","GALR2","GALR3",
    "MC1R","MC2R","MC3R","MC4R","MC5R",
    # ── Dopamine system ───────────────────────────────────────────────────────
    "TH","DDC","DBH","PNMT","DRD1","DRD2","DRD3","DRD4","DRD5",
    "SLC6A3","SLC18A1","SLC18A2","COMT","MAOA","MAOB","ALDH2",
    # ── Serotonin system ──────────────────────────────────────────────────────
    "TPH1","TPH2","SLC6A4","MAOA",
    "HTR1A","HTR1B","HTR1D","HTR1E","HTR2A","HTR2B","HTR2C",
    "HTR3A","HTR3B","HTR4","HTR5A","HTR6","HTR7",
    # ── Norepinephrine system ─────────────────────────────────────────────────
    "DBH","SLC6A2","ADRA1A","ADRA1B","ADRA2A","ADRA2B","ADRA2C",
    "ADRB1","ADRB2","ADRB3",
    # ── GABA system ───────────────────────────────────────────────────────────
    "GAD1","GAD2","GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6",
    "GABRB1","GABRB2","GABRB3","GABRG1","GABRG2","GABRG3","GABRD","GABRE",
    "GABRR1","GABRR2","GABBR1","GABBR2","SLC6A1","SLC6A11","ABAT",
    # ── Glutamate system ──────────────────────────────────────────────────────
    "GRIN1","GRIN2A","GRIN2B","GRIN2C","GRIN2D","GRIN3A","GRIN3B",
    "GRIA1","GRIA2","GRIA3","GRIA4","GRIK1","GRIK2","GRIK3","GRIK4","GRIK5",
    "GRM1","GRM2","GRM3","GRM4","GRM5","GRM6","GRM7","GRM8",
    "SLC1A1","SLC1A3","SLC1A6","SLC17A5","SLC17A7","SLC17A6",
    # ── Cholinergic system ────────────────────────────────────────────────────
    "CHAT","SLC18A3","CHRM1","CHRM2","CHRM3","CHRM4","CHRM5",
    "CHRNA7","CHRNA4","CHRNB2","CHRNA3","CHRNB4",
    # ── Opioid / cannabinoid / adenosine ──────────────────────────────────────
    "OPRM1","OPRD1","OPRK1","OPRL1","PENK","PDYN","POMC",
    "CNR1","CNR2","FAAH","MGLL","NAPEPLD",
    "ADORA1","ADORA2A","ADORA2B","ADORA3",
    # ── HPA / stress axis ─────────────────────────────────────────────────────
    "NR3C1","NR3C2","FKBP5","HSD11B1","HSD11B2","CYP11B1","CYP11B2",
    # ── Nitric oxide / NOS ────────────────────────────────────────────────────
    "NOS1","NOS2","NOS3",
    # ── Ion channels / transporters ───────────────────────────────────────────
    "SCN1A","SCN2A","SCN8A","KCNQ2","KCNQ3","KCNJ10","KCNJ6",
    "HCN1","HCN2","HCN4","CACNA1A","CACNA1C","CACNA1D","CACNA2D1",
    # ── Circadian ────────────────────────────────────────────────────────────
    "CLOCK","PER1","PER2","PER3","CRY1","CRY2","ARNTL","ARNTL2","NR1D1","NR1D2","TIMELESS",
    # ── Transcription factors ─────────────────────────────────────────────────
    "CREB1","CREB3L2","ATF3","ATF4","ATF6",
    "FOS","FOSB","FOSL1","FOSL2","JUN","JUNB","JUND",
    "EGR1","EGR2","EGR3","EGR4","NPAS4","NPAS2","MEF2C","MEF2A","MEF2D",
    "NRF1","NFE2L2","SP1","TP53","CEBPA","CEBPB","CEBPD",
    # ── Synaptic plasticity / scaffolding ────────────────────────────────────
    "SYN1","SYN2","SYN3","SYT1","SYT7","VAMP1","VAMP2",
    "DLG1","DLG2","DLG3","DLG4","SHANK1","SHANK2","SHANK3",
    "HOMER1","HOMER2","HOMER3","GRIP1","SYNGAP1","DLGAP1","DLGAP2",
    "NCAM1","NCAM2","NRXN1","NRXN2","NRXN3",
    "NLGN1","NLGN2","NLGN3","NLGN4X",
    # ── Neuroinflammation / cytokines ─────────────────────────────────────────
    "TNF","IL1B","IL2","IL4","IL5","IL6","IL7","IL10","IL12A","IL12B",
    "IL13","IL15","IL17A","IL18","IL21","IL23A","IL33",
    "IFNG","IFNA1","IFNB1","TGFB1","TGFB2","TGFB3",
    "CXCL8","CXCL1","CXCL2","CXCL9","CXCL10","CXCL12","CXCL13",
    "CCL2","CCL3","CCL4","CCL5","CCL11","CCL20",
    "TNFSF10","TNFRSF1A","TNFRSF1B",
    "PTGS1","PTGS2","NOS2","MPO","HMGB1",
    # ── JAK-STAT / NF-kB signaling ────────────────────────────────────────────
    "STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6",
    "JAK1","JAK2","JAK3","TYK2",
    "NFKB1","NFKB2","RELA","RELB","IKBKB","IKBKG",
    "IRF1","IRF3","IRF7","MYD88",
    # ── TLR / innate immunity ─────────────────────────────────────────────────
    "TLR1","TLR2","TLR3","TLR4","TLR5","TLR7","TLR9",
    "NLRP1","NLRP3","NLRC4","CASP1","PYCARD",
    "NOD1","NOD2",
    # ── MAPK / PI3K / mTOR signaling ─────────────────────────────────────────
    "MAPK1","MAPK3","MAPK8","MAPK9","MAPK14",
    "PIK3CA","PIK3CB","PIK3CD","AKT1","AKT2","AKT3","MTOR","PTEN","TSC1","TSC2",
    "BRAF","RAF1","MAP2K1","MAP2K2",
    "GSK3A","GSK3B","CAMK2A","CAMK2B","CAMK4","CDK5","DYRK1A",
    # ── Apoptosis / cell survival ─────────────────────────────────────────────
    "BCL2","BAX","BCL2L1","MCL1","BID","BAD","PMAIP1",
    "CASP3","CASP8","CASP9","CASP1","APAF1","CYCS",
    "TP53","CDKN1A","CDKN2A","MDM2",
    # ── Oxidative stress / antioxidant ───────────────────────────────────────
    "SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","PRDX2","TXN","TXNRD1",
    "HIF1A","EPAS1","NFE2L2","HMOX1","HSPA5","HSP90AA1","HSPB1",
    # ── Autophagy / ubiquitin ─────────────────────────────────────────────────
    "BECN1","ATG5","ATG7","ATG12","ATG16L1","SQSTM1","OPTN","CALCOCO2",
    "PINK1","PRKN","LRRK2","ULK1","LAMP2",
    # ── DNA damage / repair / epigenetic ─────────────────────────────────────
    "DNMT1","DNMT3A","DNMT3B","TET1","TET2","TET3",
    "HDAC1","HDAC2","HDAC3","HDAC5","HDAC6","SIRT1","SIRT2","SIRT3",
    "KAT2A","KAT2B","KAT6A","EP300","CREBBP",
    "EZH2","SUZ12","EED","RING1","BMI1","KDM6A","KDM5C",
    "MECP2","RETT","MBD1","MBD2","MBD3",
    # ── Glial / astrocyte / oligodendrocyte ──────────────────────────────────
    "GFAP","S100B","AQP4","ALDH1L1","GLUL","GJA1","GJB6","VIM",
    "MOG","MBP","PLP1","CNP","MAG","MOBP","MYRF","OLIG1","OLIG2",
    "NKX2-2","SOX10",
    # ── Microglia ─────────────────────────────────────────────────────────────
    "AIF1","CX3CR1","CX3CL1","TMEM119","P2RY12","P2RY13","SALL1","HEXB",
    "CD68","TREM2","ITGAM","PTPRC","CSF1R","SIGLECH",
    # ── Alzheimer's-specific ──────────────────────────────────────────────────
    "APP","PSEN1","PSEN2","APOE","MAPT","CLU","BIN1","CR1","ABCA7",
    "BACE1","BACE2","ADAM10","ADAM17","APOJ",
    "TREM2","CLU","PICALM","CD33","MS4A6A","EPHA1","CD2AP","FERMT2",
    # ── Parkinson's-specific ──────────────────────────────────────────────────
    "SNCA","LRRK2","PRKN","PINK1","PARK7","GBA","UCHL1","ATP13A2",
    "FBXO7","PLA2G6","DNAJC6","SYNJ1","VPS35","EIF4G1","CHCHD2",
    # ── Autism / ASD ──────────────────────────────────────────────────────────
    "SHANK1","SHANK2","SHANK3","SYNGAP1","NRXN1","NRXN2","NRXN3",
    "NLGN1","NLGN2","NLGN3","NLGN4X","CNTNAP2","CNTNAP4",
    "CHD8","DYRK1A","SCN2A","SCN1A","ADNP","PTEN","TSC1","TSC2",
    "FMR1","MECP2","RETT","GABRB3","FOXP1","FOXP2","TCF4","ANK2","ANK3",
    # ── Schizophrenia-specific ────────────────────────────────────────────────
    "DISC1","NRG1","DTNBP1","RELN","DAOA","G72","COMT","PRODH",
    "NRXN1","CNTNAP2","DISC2",
    # ── ADHD-specific ─────────────────────────────────────────────────────────
    "LPHN3","SNAP25","DRD4","DRD5","SLC6A4","SLC6A3","HTR1B","BDNF","NOS1",
    # ── PTSD-specific ─────────────────────────────────────────────────────────
    "FKBP5","NR3C2","ADCYAP1R1","CRHR1","CRHR2","AVP","OXTR",
    # ── OCD-specific ──────────────────────────────────────────────────────────
    "SLC6A4","HTR2A","HTR2C","DRD4","SAPAP3","DLGAP3",
    # ── General signaling / Misc ──────────────────────────────────────────────
    "ACE","ACE2","BDKRB1","BDKRB2","LEPR","INSR","PPARA","PPARG","RXRA",
    "NR3C1","NR3C2","ESR1","ESR2","AR","THRB","RARA",
    "NCAM1","NCAM2","L1CAM","NRCAM","CDON","BOC",
    "SLC6A4","SLC6A3","SLC6A2","SLC6A1",
    "IDO1","IDO2","TDO2","KYNU","HAAO",
    "HMGB1","CRP","SAA1","SAA2",
    "RET","PHOX2B","PHOX2A",
    "TJP1","FFAR2","VIPR1","GLP1R","GLP2R",
    "CHD8","ANK2","TCF4","CNTNAP2","FOXP1","DYRK1A","SCN2A",
    "CAMK2A","NRGN","SYN1","NCAM1","SNCA","BECN1","SQSTM1","ACE2",
]

GUT_GENES = [
    # ── Intestinal stem cell / epithelium ─────────────────────────────────────
    "VIL1","CDX2","LGR5","MUC2","MUC1","MUC5AC","MUC5B","MUC4","MUC13",
    "ATOH1","DEFA5","DEFA6","CHGA","CHGB","GIP",
    "SLC5A1","FABP1","FABP2","FABP6","SI","EPCAM","OLFM4","REG3A","REG1A",
    "AQP3","AQP8","SLC26A3","SLC26A2","GLP1R","PYY",
    "ASCL2","SMOC2","LYZ","NEUROG3","PAX4","INSM1","GCG","SCT",
    "GHRL","TPH1","CA2","CDH1","TFF1","TFF2","TFF3",
    "FCGBP","ANPEP","SLC15A1","APOA4","APOB","ALPI","SLC2A5","LEP",
    "SLC10A2","ABCG2","ABCB11","TRPV1","TRPV4","TRPV3",
    # ── Tight junction / intestinal barrier ───────────────────────────────────
    "OCLN","TJP1","TJP2","TJP3",
    "CLDN1","CLDN2","CLDN3","CLDN4","CLDN5","CLDN6","CLDN7","CLDN8",
    "CLDN10","CLDN11","CLDN12","CLDN14","CLDN15","CLDN18","CLDN19",
    "JAM1","F11R","JAM2","JAM3","MARVELD2","MARVELD3","MARVELD1",
    "AFDN","MLLT4","ESAM","INADL",
    # ── Gut inflammation / immune ──────────────────────────────────────────────
    "IL1B","IL1A","IL2","IL4","IL5","IL6","IL8","IL10","IL12A","IL12B",
    "IL13","IL17A","IL17F","IL18","IL22","IL23A","IL25","IL33","IL34",
    "IFNG","IFNA1","IFNB1","IFNL1","TNF","TGFB1","TGFB2","TGFB3",
    "CXCL8","CXCL1","CXCL2","CXCL9","CXCL10","CXCL12",
    "CCL2","CCL3","CCL4","CCL5","CCL11","CCL20","CCL25","CCL28",
    "PTGS2","PTGS1","NOS2","MPO","CALPR","S100A8","S100A9","S100A12",
    "TLR2","TLR3","TLR4","TLR5","TLR7","TLR9",
    "NOD1","NOD2","NLRP3","NLRP6","PYCARD","CASP1","IL1RN",
    "MUC2","FCGR2A","FCGR3A","CD14","LBP","HMGB1",
    "TNFSF10","TNFRSF1A","TNFRSF1B","TNFSF13B",
    "STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B",
    "NFKB1","NFKB2","RELA","IKBKB","IRF1","IRF3","MYD88",
    # ── Enteric nervous system / serotonin ────────────────────────────────────
    "TPH1","SLC6A4","HTR3A","HTR4","HTR1A","HTR2B","HTR7","CHRM3","NOS1","NOS3",
    "VIP","VIPR1","VIPR2","NPY","PENK","TAC1","GAL","ADCYAP1","CALCA",
    "GDNF","NTF3","BDNF","NGF","NTRK2","RET","PHOX2B","PHOX2A",
    "ADORA1","ADORA2A","P2RX3","P2RY1","P2RY4",
    # ── Gut hormones / signaling ──────────────────────────────────────────────
    "NR3C1","CRH","LEPR","INSR","CCK","AGRP","FFAR2","FFAR3","FFAR4",
    "DRD2","DRD4","TLR4","HMGB1","CHAT","IDO1","IDO2","TDO2","KYNU",
    "CLOCK","PER2","ARNTL","NR1D1","CRY1","CRY2",
    "GLP1R","GLP2R","GIPR","GCGR","SSTR1","SSTR2","SSTR3","SSTR5",
    "PPARG","PPARA","RXRA","HIF1A","NFE2L2","HMOX1",
    # ── Parkinson's / neurodegeneration in gut ────────────────────────────────
    "SNCA","LRRK2","PRKN","PINK1","BECN1","SQSTM1","ACE2",
    # ── Autism / ASD gut-expressed genes ─────────────────────────────────────
    "SHANK3","CHD8","CNTNAP2","NRXN1","NLGN3","SYNGAP1","TCF4","ANK2",
    "SYN1","NCAM1","FOXP1","FMR1","MECP2","TSC1","TSC2","PTEN",
    # ── Oxidative stress / barrier integrity ──────────────────────────────────
    "SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","TXNRD1",
    "NFE2L2","HMOX1","HSP90AA1","HSPA5","HSPB1",
    "ACSL4","LPCAT3","GPX4",
    # ── Gut microbiome interaction ────────────────────────────────────────────
    "REG1A","REG3A","DEFB1","DEFB4A","DEFA5","DEFA6","LYZ",
    "PIGR","JCHAIN","IGHM","IGHA1","IGHA2","IL15","TSLP",
    # ── Motility / muscle ─────────────────────────────────────────────────────
    "ACTG2","MYH11","MYLK","CNN1","DES","VIM",
    # ── Mucus / secretion ─────────────────────────────────────────────────────
    "MUC2","MUC5B","MUC5AC","TFF1","TFF3","CLCA1","SPDEF",
    # ── Absorption / metabolism ────────────────────────────────────────────────
    "SLC6A19","SLC7A9","SLC6A14","SLC16A1","ABCB1","ABCC2","ABCC3",
    "CYP3A4","CYP3A5","CYP2C9","UGT1A1","SULT1A1",
]

print("Dynamic gene detection mode: any gene symbol found in biological context will be captured.")


# ─────────────────────────────────────────────────────────────────────────────
# PARSING AND DISAMBIGUATION SUPPORT
# ─────────────────────────────────────────────────────────────────────────────

_BIO_CONTEXT_RE = re.compile(
    r"\b(gene|genes|expression|expressed|mrna|protein|receptor|knockout|"
    r"knockdown|overexpression|transcript|transcription|encoded|encodes|"
    r"regulated|upregulated|downregulated|levels|signaling|signalling|"
    r"pathway|mutation|variant|allele|promoter|locus|immunoreactivity|"
    r"immunostaining|antibody|antibodies|staining|positive|neurons|"
    r"neuronal|intestinal|colonic|hippocampal|cortical|striatal|"
    r"deficiency|deficient|null mutant|heterozygous|homozygous|transgenic|"
    r"knockin|overexpressing|silencing|inhibition|activation)\b", re.IGNORECASE
)

_TIER2_ALIASES = {
    # Gut disambiguation
    "SI":     re.compile(r"\b(sucrase.isomaltase|sucrase|isomaltase|disaccharidase|brush.border enzyme|brush border membrane)\b", re.IGNORECASE),
    "GIP":    re.compile(r"\b(gastric inhibitory polypeptide|glucose.dependent insulinotropic|incretin|gip receptor|gipr|k.cell|gip secretion)\b", re.IGNORECASE),
    # Tight junction / barrier — common names
    "OCLN":   re.compile(r"\boccludin\b", re.IGNORECASE),
    "TJP1":   re.compile(r"\b(ZO[\s\-]?1|zona occludens[\s\-]?1|zonula occludens[\s\-]?1)\b", re.IGNORECASE),
    "TJP2":   re.compile(r"\b(ZO[\s\-]?2|zona occludens[\s\-]?2|zonula occludens[\s\-]?2)\b", re.IGNORECASE),
    "CLDN1":  re.compile(r"\bclaudin[\s\-]?1\b", re.IGNORECASE),
    "CLDN2":  re.compile(r"\bclaudin[\s\-]?2\b", re.IGNORECASE),
    "CLDN3":  re.compile(r"\bclaudin[\s\-]?3\b", re.IGNORECASE),
    "CLDN4":  re.compile(r"\bclaudin[\s\-]?4\b", re.IGNORECASE),
    "CLDN5":  re.compile(r"\bclaudin[\s\-]?5\b", re.IGNORECASE),
    "CLDN7":  re.compile(r"\bclaudin[\s\-]?7\b", re.IGNORECASE),
    "CLDN8":  re.compile(r"\bclaudin[\s\-]?8\b", re.IGNORECASE),
    "CLDN12": re.compile(r"\bclaudin[\s\-]?12\b", re.IGNORECASE),
    "CLDN15": re.compile(r"\bclaudin[\s\-]?15\b", re.IGNORECASE),
    # Inflammatory cytokines — common names
    "IL1B":   re.compile(r"\b(IL[\s\-]?1[\s\-]?(?:beta|β|b)\b|interleukin[\s\-]?1[\s\-]?(?:beta|β|b))\b", re.IGNORECASE),
    "IL6":    re.compile(r"\b(IL[\s\-]?6\b|interleukin[\s\-]?6)\b", re.IGNORECASE),
    "IL10":   re.compile(r"\b(IL[\s\-]?10\b|interleukin[\s\-]?10)\b", re.IGNORECASE),
    "IL17A":  re.compile(r"\b(IL[\s\-]?17[A]?\b|interleukin[\s\-]?17)\b", re.IGNORECASE),
    "IL18":   re.compile(r"\b(IL[\s\-]?18\b|interleukin[\s\-]?18)\b", re.IGNORECASE),
    "IL22":   re.compile(r"\b(IL[\s\-]?22\b|interleukin[\s\-]?22)\b", re.IGNORECASE),
    "IFNG":   re.compile(r"\b(IFN[\s\-]?(?:gamma|γ|g)\b|interferon[\s\-]?(?:gamma|γ))\b", re.IGNORECASE),
    "TNF":    re.compile(r"\b(TNF[\s\-]?(?:alpha|α|a)\b|tumor\s+necrosis\s+factor[\s\-]?(?:alpha|α)?)\b", re.IGNORECASE),
    "TGFB1":  re.compile(r"\b(TGF[\s\-]?(?:beta|β|b)[\s\-]?1?\b|transforming\s+growth\s+factor[\s\-]?(?:beta|β))\b", re.IGNORECASE),
    "CXCL8":  re.compile(r"\b(IL[\s\-]?8\b|interleukin[\s\-]?8|CXCL[\s\-]?8)\b", re.IGNORECASE),
    "CCL2":   re.compile(r"\b(MCP[\s\-]?1\b|monocyte\s+chemoattractant\s+protein[\s\-]?1|CCL[\s\-]?2)\b", re.IGNORECASE),
    # Inflammasome
    "NLRP3":  re.compile(r"\b(NLRP[\s\-]?3|NALP[\s\-]?3|cryopyrin)\b", re.IGNORECASE),
    # Enzyme common names
    "PTGS2":  re.compile(r"\b(COX[\s\-]?2|cyclooxygenase[\s\-]?2)\b", re.IGNORECASE),
    "NOS2":   re.compile(r"\b(iNOS|inducible\s+NOS|inducible\s+nitric\s+oxide\s+synthase)\b", re.IGNORECASE),
    # Mucins
    "MUC2":   re.compile(r"\b(mucin[\s\-]?2|MUC[\s\-]?2)\b", re.IGNORECASE),
    "MUC1":   re.compile(r"\b(mucin[\s\-]?1|MUC[\s\-]?1)\b", re.IGNORECASE),
    # Brain genes common names
    "BDNF":   re.compile(r"\b(brain[\s\-]?derived\s+neurotrophic\s+factor|BDNF)\b", re.IGNORECASE),
    "GDNF":   re.compile(r"\b(glial\s+cell[\s\-]?(?:line[\s\-]?derived)?\s+neurotrophic\s+factor|GDNF)\b", re.IGNORECASE),
    "SNCA":   re.compile(r"\b(alpha[\s\-]?synuclein|α[\s\-]?synuclein|SNCA)\b", re.IGNORECASE),
    "MAPT":   re.compile(r"\b(tau\s+protein|hyperphosphorylated\s+tau|p[\s\-]?tau|tau\s+pathology)\b", re.IGNORECASE),
    "APP":    re.compile(r"\b(amyloid\s+precursor\s+protein|amyloid[\s\-]?β|amyloid[\s\-]?beta|Aβ\s*\d*|A-beta|beta[\s\-]?amyloid)\b", re.IGNORECASE),
    "APOE":   re.compile(r"\b(apolipoprotein\s*E|ApoE[\s\-]?\d?|APOE[\s\-]?\d?)\b", re.IGNORECASE),
    "FKBP5":  re.compile(r"\b(FKBP51|FK506[\s\-]?binding\s+protein\s+5)\b", re.IGNORECASE),
    "NFE2L2": re.compile(r"\b(NRF[\s\-]?2|nuclear\s+factor\s+erythroid\s+2|NF[\s\-]?E2L2)\b", re.IGNORECASE),
    "PPARGC1A": re.compile(r"\b(PGC[\s\-]?1[\s\-]?(?:alpha|α)?|peroxisome\s+proliferator[\s\-]?activated\s+receptor\s+gamma\s+coactivator)\b", re.IGNORECASE),
    "HMOX1":  re.compile(r"\b(HO[\s\-]?1|heme\s+oxygenase[\s\-]?1)\b", re.IGNORECASE),
    # Additional cytokine aliases
    "IL2":    re.compile(r"\b(IL[\s\-]?2\b|interleukin[\s\-]?2)\b", re.IGNORECASE),
    "IL4":    re.compile(r"\b(IL[\s\-]?4\b|interleukin[\s\-]?4)\b", re.IGNORECASE),
    "IL12A":  re.compile(r"\b(IL[\s\-]?12[A]?\b|interleukin[\s\-]?12)\b", re.IGNORECASE),
    "IL13":   re.compile(r"\b(IL[\s\-]?13\b|interleukin[\s\-]?13)\b", re.IGNORECASE),
    "IL33":   re.compile(r"\b(IL[\s\-]?33\b|interleukin[\s\-]?33)\b", re.IGNORECASE),
    "CXCL10": re.compile(r"\b(IP[\s\-]?10\b|CXCL[\s\-]?10|interferon[\s\-]?gamma[\s\-]?induced\s+protein\s+10)\b", re.IGNORECASE),
    "CXCL12": re.compile(r"\b(SDF[\s\-]?1|stromal\s+cell[\s\-]?derived\s+factor[\s\-]?1|CXCL[\s\-]?12)\b", re.IGNORECASE),
    "CCL5":   re.compile(r"\b(RANTES\b|CCL[\s\-]?5)\b", re.IGNORECASE),
    # Receptor aliases
    "NFKB1":  re.compile(r"\b(NF[\s\-]?κ?B[\s\-]?p50|NF[\s\-]?kB|nuclear\s+factor[\s\-]?kappa\s*B)\b", re.IGNORECASE),
    "ADCYAP1R1": re.compile(r"\b(PAC1\s*R?|PAC[\s\-]?1\s+receptor|PACAP\s+receptor|PAC1[\s\-]?receptor)\b", re.IGNORECASE),
    "OPRM1":  re.compile(r"\b(mu[\s\-]?opioid\s+receptor|MOR[\s\-]?1|MOP\s+receptor)\b", re.IGNORECASE),
    "CNR1":   re.compile(r"\b(CB[\s\-]?1\s+receptor|cannabinoid\s+receptor\s+1|CB1R)\b", re.IGNORECASE),
    "CNR2":   re.compile(r"\b(CB[\s\-]?2\s+receptor|cannabinoid\s+receptor\s+2|CB2R)\b", re.IGNORECASE),
    "FAAH":   re.compile(r"\b(fatty\s+acid\s+amide\s+hydrolase|FAAH)\b", re.IGNORECASE),
    "DISC1":  re.compile(r"\b(disrupted\s+in\s+schizophrenia[\s\-]?1|DISC[\s\-]?1)\b", re.IGNORECASE),
    "NRG1":   re.compile(r"\b(neuregulin[\s\-]?1|heregulin|NRG[\s\-]?1)\b", re.IGNORECASE),
    "RELN":   re.compile(r"\b(reelin|RELN)\b", re.IGNORECASE),
    "MAOA":   re.compile(r"\b(monoamine\s+oxidase[\s\-]?A|MAO[\s\-]?A)\b", re.IGNORECASE),
    "MAOB":   re.compile(r"\b(monoamine\s+oxidase[\s\-]?B|MAO[\s\-]?B)\b", re.IGNORECASE),
    "COMT":   re.compile(r"\b(catechol[\s\-]?O[\s\-]?methyltransferase|COMT)\b", re.IGNORECASE),
    "GBA":    re.compile(r"\b(glucocerebrosidase|beta[\s\-]?glucocerebrosidase|GCase)\b", re.IGNORECASE),
    "FMR1":   re.compile(r"\b(fragile\s+X\s+(?:mental\s+retardation\s+protein|syndrome\s+protein|protein)|FMRP)\b", re.IGNORECASE),
    "MECP2":  re.compile(r"\b(methyl[\s\-]?CpG[\s\-]?binding\s+protein\s+2|MeCP2|MECP[\s\-]?2)\b", re.IGNORECASE),
    "TJP3":   re.compile(r"\b(ZO[\s\-]?3|zona\s+occludens[\s\-]?3)\b", re.IGNORECASE),
    # Claudin extensions
    "CLDN6":  re.compile(r"\bclaudin[\s\-]?6\b", re.IGNORECASE),
    "CLDN10": re.compile(r"\bclaudin[\s\-]?10\b", re.IGNORECASE),
    "CLDN11": re.compile(r"\bclaudin[\s\-]?11\b", re.IGNORECASE),
    "CLDN14": re.compile(r"\bclaudin[\s\-]?14\b", re.IGNORECASE),
    "CLDN18": re.compile(r"\bclaudin[\s\-]?18\b", re.IGNORECASE),
    "CLDN19": re.compile(r"\bclaudin[\s\-]?19\b", re.IGNORECASE),
    # Gut barrier / permeability proteins by common name
    "TJP1":   re.compile(r"\b(ZO[\s\-]?1|zona\s+occludens[\s\-]?1|tight[\s\-]?junction\s+protein[\s\-]?1)\b", re.IGNORECASE),
    "MUC2":   re.compile(r"\b(mucin[\s\-]?2|goblet[\s\-]?cell\s+mucin)\b", re.IGNORECASE),
    "MUC5AC": re.compile(r"\bmucin[\s\-]?5AC\b", re.IGNORECASE),
    "FFAR2":  re.compile(r"\b(free\s+fatty\s+acid\s+receptor[\s\-]?2|GPR43|SCFA\s+receptor)\b", re.IGNORECASE),
    "FFAR3":  re.compile(r"\b(free\s+fatty\s+acid\s+receptor[\s\-]?3|GPR41)\b", re.IGNORECASE),
    "GLP1R":  re.compile(r"\b(GLP[\s\-]?1\s+receptor|glucagon[\s\-]?like\s+peptide[\s\-]?1\s+receptor)\b", re.IGNORECASE),
    "GCG":    re.compile(r"\b(glucagon[\s\-]?like\s+peptide[\s\-]?1\b|GLP[\s\-]?1\b|proglucagon)\b", re.IGNORECASE),
    "GHRL":   re.compile(r"\b(ghrelin|acyl[\s\-]?ghrelin|des[\s\-]?acyl[\s\-]?ghrelin)\b", re.IGNORECASE),
    "PYY":    re.compile(r"\b(peptide\s+YY|peptide[\s\-]?tyrosine[\s\-]?tyrosine)\b", re.IGNORECASE),
    "CCK":    re.compile(r"\b(cholecystokinin|CCK[\s\-]?\d*)\b", re.IGNORECASE),
    "VIL1":   re.compile(r"\bvillin[\s\-]?1?\b", re.IGNORECASE),
    "CHGA":   re.compile(r"\b(chromogranin[\s\-]?A|CgA)\b", re.IGNORECASE),
    "TPH1":   re.compile(r"\b(tryptophan\s+hydroxylase[\s\-]?1|TPH[\s\-]?1)\b", re.IGNORECASE),
    "IDO1":   re.compile(r"\b(indoleamine[\s\-]?2,?3[\s\-]?dioxygenase[\s\-]?1?|IDO[\s\-]?1)\b", re.IGNORECASE),
    "KYNU":   re.compile(r"\b(kynureninase|kynurenine\s+aminotransferase)\b", re.IGNORECASE),
    "ACE2":   re.compile(r"\b(angiotensin[\s\-]?converting\s+enzyme[\s\-]?2|ACE[\s\-]?2)\b", re.IGNORECASE),
    "DEFB1":  re.compile(r"\b(defensin[\s\-]?beta[\s\-]?1|beta[\s\-]?defensin[\s\-]?1|hBD[\s\-]?1)\b", re.IGNORECASE),
    "LYZ":    re.compile(r"\b(lysozyme|muramidase)\b", re.IGNORECASE),
}

BANNED_ACRONYM_CONTEXT = {
    "SI": [r"suicidal", r"ideation", r"severity", r"international system", r"signal intensity", r"saturation index", r"substantia innominata"],
    "TH": [r"\bth\s+cell", r"\bth\d+\b", r"\bth\s+hormone", r"t\s*cell helper", r"helper\s*t", r"thyroid", r"tuesday", r"thursday", r"\b\d+th\b"]
}

_NEUROPEPTIDE_INTESTINAL_RE = re.compile(
    r"\bvasoactive\s+intestinal\s+(?:peptide|polypeptide)\b|"
    r"\bVIP[\s\-](?:expressing|positive|\+|ir|immunoreactive|ergic|neurons?|interneurons?|cells?)\b|"
    r"\bVIP\+\b",
    re.IGNORECASE,
)

_BRAIN_KW_RE = re.compile(r"\b(?:" + "|".join(re.escape(kw) for kw in BRAIN_KEYWORDS) + r")\b", re.IGNORECASE)
_GUT_KW_RE   = re.compile(r"\b(?:" + "|".join(re.escape(kw) for kw in GUT_KEYWORDS)   + r")\b", re.IGNORECASE)

def _mask_neuropeptide_names(text):
    """Remove intestinal-containing neuropeptide names so they don't trigger gut classification."""
    return _NEUROPEPTIDE_INTESTINAL_RE.sub("NEUROPEPTIDE_NAME", text)

def detect_tissue(text):
    masked = _mask_neuropeptide_names(text)
    brain_found = bool(_BRAIN_KW_RE.search(masked))
    gut_found   = bool(_GUT_KW_RE.search(masked))
    return "Brain & Gut" if (brain_found and gut_found) else ("Brain" if brain_found else ("Gut" if gut_found else "Other"))

def detect_tissue_for_gene(gene, text):
    """Determine tissue context for a specific gene based on WHERE it is studied, not just mentioned.

    Strategy: check results section first. If the gene appears there with a tissue context,
    that is the definitive answer — background/intro mentions of other tissues are ignored.
    Only fall back to full-abstract scanning if the gene doesn't appear in the results section.
    """
    results_text = _get_results_section(text)
    sentences    = re.split(r'(?<=[.!?])\s+', text)
    results_sentences = re.split(r'(?<=[.!?])\s+', results_text)
    gene_pat = re.compile(r'\b' + re.escape(gene) + r'\b', re.IGNORECASE)

    _FLUID_KW_RE = re.compile(
        r"\b(?:serum|plasma|blood|CSF|cerebrospinal\s+fluid|urine|urinary|saliva|salivary|"
        r"peripheral\s+blood|whole\s+blood|venous\s+blood)\b", re.IGNORECASE,
    )
    # Sentences describing bacteria/microbiota should not trigger gut tissue classification
    _BACTERIA_KW_RE = re.compile(
        r"\b(?:bacteria(?:l)?|microbiota|microbiome|microorganism|16[Ss]\s+rRNA|"
        r"fecal\s+microbiota|gut\s+microbiota|gut\s+microbiome|microbial\s+(?:composition|community|"
        r"diversity|abundance|profile)|Lactobacillus|Bifidobacterium|Firmicutes|Bacteroidetes|"
        r"Akkermansia|Faecalibacterium|Clostridium|Prevotella|Ruminococcus)\b", re.IGNORECASE,
    )

    def _tissue_in_sentences(sent_list):
        brain, gut = False, False
        for i, sent in enumerate(sent_list):
            if not gene_pat.search(sent):
                continue
            # Skip sentences where the gene is only mentioned in fluid context
            if _FLUID_KW_RE.search(sent) and not _BRAIN_KW_RE.search(sent) and not _GUT_KW_RE.search(sent):
                continue
            # Skip sentences where "gut" context is actually gut bacteria, not gut tissue
            if _BACTERIA_KW_RE.search(sent) and not _BRAIN_KW_RE.search(sent):
                continue
            masked_sent = _mask_neuropeptide_names(sent)
            sent_brain = bool(_BRAIN_KW_RE.search(masked_sent))
            sent_gut   = bool(_GUT_KW_RE.search(masked_sent))
            if sent_brain or sent_gut:
                # The gene's own sentence has a tissue signal — use it directly
                if sent_brain: brain = True
                if sent_gut:   gut   = True
            else:
                # No tissue in the gene sentence — expand to ±1 sentence only
                window_start = max(0, i - 1)
                window_end   = min(len(sent_list), i + 2)
                context = " ".join(sent_list[window_start:window_end])
                masked  = _mask_neuropeptide_names(context)
                if _BRAIN_KW_RE.search(masked): brain = True
                if _GUT_KW_RE.search(masked):   gut   = True
        return brain, gut

    # 1. Check results section only
    brain_res, gut_res = _tissue_in_sentences(results_sentences)
    if brain_res or gut_res:
        return ("Brain & Gut" if (brain_res and gut_res)
                else ("Brain" if brain_res else "Gut"))

    # 2. Fall back to full abstract (gene only mentioned in background)
    brain_all, gut_all = _tissue_in_sentences(sentences)
    if brain_all or gut_all:
        return ("Brain & Gut" if (brain_all and gut_all)
                else ("Brain" if brain_all else "Gut"))

    # 3. Last resort: abstract-level tissue
    return detect_tissue(text)

def detect_species(text):
    tl = text.lower()
    first_half = tl[:len(tl)//2]
    found, rodents = set(), set()
    for species, kws in RODENT_SPECIES.items():
        for kw in kws:
            if re.search(r"\b" + re.escape(kw) + r"\b", tl):
                found.add(species); rodents.add(species); break
    if re.search(r"\brodents?\b", tl) and not rodents:
        found.add("Other Rodents")
    for species, kws in SPECIES_KEYWORDS.items():
        st = first_half if species == "Humans" else tl
        for kw in kws:
            if re.search(r"\b" + re.escape(kw) + r"\b", st):
                found.add(species); break
    return list(found)

def detect_gene_expr(text):    return any(p.search(text) for p in GENE_EXPR_PATTERNS)
def detect_microbiome(text):   return any(p.search(text) for p in MICROBIOME_PATTERNS)
_LINKING_WORDS = {
    "of", "in", "was", "is", "were", "are", "the", "its", "their",
    "a", "an", "by", "for", "that", "which", "to", "be", "been",
}

_UP_RE   = re.compile(
    r"\b(?:upregulat\w+|up[\s\-]regulat\w+|over[\s\-]?express\w+|"
    r"increas\w+|elevat\w+|higher|enhanc\w+|augment\w+|amplif\w+|"
    r"potentiat\w+|activat\w+|induced\b|induction\b|"
    r"gain[\s\-]of[\s\-]function|positively\s+regulat\w+|"
    r"\d+[\s\-]?fold\s+(?:increase|higher|elevation)|"
    r"greater\s+(?:expression|levels?)|excess\w+)\b",
    re.IGNORECASE
)
_DOWN_RE = re.compile(
    r"\b(?:downregulat\w+|down[\s\-]regulat\w+|under[\s\-]?express\w+|"
    r"decreas\w+|reduc\w+|lower|suppress\w+|attenuate\w+|"
    r"diminish\w+|blunt\w+|abolish\w+|deplet\w+|"
    r"loss\s+of\s+expression|impair\w+\s+expression|"
    r"inhibit\w+|repress\w+|silenc\w+|knock(?:ed|ing)?\s*down|knockdown|"
    r"negatively\s+regulat\w+|deficien\w+\s+(?:in|of)\s+\w+\s+expression|"
    r"\d+[\s\-]?fold\s+(?:decrease|lower|reduction)|"
    r"absent\s+expression|loss\s+of\s+\w+\s+(?:expression|protein|mRNA))\b",
    re.IGNORECASE
)

def detect_gene_regulation(gene, text):
    """Return 'Upregulated', 'Downregulated', 'Mixed', or 'Unknown' for a gene
    based on regulation language in the same or adjacent sentence."""
    results_text = _get_results_section(text)
    if not results_text or len(results_text) < 50:
        results_text = text
    sentences    = re.split(r'(?<=[.!?])\s+', results_text)
    gene_pat     = re.compile(r'\b' + re.escape(gene) + r'\b', re.IGNORECASE)
    # also check alias
    alias_pat = _TIER2_ALIASES.get(gene)

    up, down = 0, 0
    for i, sent in enumerate(sentences):
        hit = gene_pat.search(sent) or (alias_pat and alias_pat.search(sent))
        if not hit:
            continue
        window = " ".join(sentences[max(0,i-1):i+2])
        if _UP_RE.search(window):   up   += 1
        if _DOWN_RE.search(window): down += 1

    if up > 0 and down > 0: return "Mixed"
    if up > 0:               return "Upregulated"
    if down > 0:             return "Downregulated"
    return "No directional language"

def _regex_prefilter(text):
    """Cheap regex pre-filter for the LLM qualification path (Option B).

    True iff the abstract is NOT excluded by a confound filter AND has a real
    disease model AND an expression-measurement signal in the results section.
    This is exactly parts 1+2 of detect_gene_context; the over-strict 1-3-word
    proximity heuristic (part 3) is deliberately NOT applied here, so the LLM
    can rescue papers whose expression signal is phrased without tight
    gene-term adjacency (the regex's main false-negative source). The cheap
    bulk-excluders stay so the LLM set stays bounded (~hours, not days).
    """
    if (
        _abstract_uses_therapeutic_drug(text)
        or _abstract_uses_nondrug_intervention(text)
        or _abstract_attributes_expression_to_treatment(text)
        or _abstract_uses_epigenetics(text)
        or _abstract_is_behavioral_primary(text)
        or _is_genetic_association_study(text)
        or _is_computational_study(text)
        or _is_fluid_only_measurement(text)
        or _is_blood_cell_line_study(text)
    ):
        return False
    if not _has_disease_model(text):
        return False
    if not _EXPRESSION_MEASUREMENT_RE.search(_get_results_section(text)):
        return False
    return True


def detect_gene_context(text, mentioned_genes):
    """Return True when a regulation/expression term is within 1–3 words of
    'gene'/'genes' or a detected gene name, measured from the *actual* regex
    match position (via finditer) so neighbouring tokens cannot accidentally
    trigger the check.

    Proximity rules (same for 'gene'/'genes' and specific gene names):
      1 word  : direct adjacency        — 'BDNF expression', 'gene expression'
      2 words : one linking word        — 'expression of BDNF', 'BDNF was upregulated'
      3 words : at least one linking    — 'SNCA and BDNF were both expressed'

    Examples that pass:  'BDNF expression', 'expression of BDNF',
                         'gene expression', 'BDNF was upregulated'
    Examples that fail:  'protein expression', 'respective protein expression',
                         'genes have been associated with their respective protein expression'
    """
    # Abstract-level exclusion + positive requirements (parts 1+2). Shared with
    # _regex_prefilter so the regex fallback and the LLM pre-filter never diverge.
    if not _regex_prefilter(text):
        return False

    # Combine 'gene'/'genes', detected gene symbols, AND their common-name aliases
    gene_terms = [r"genes?"] + [re.escape(g) for g in mentioned_genes]
    for g in mentioned_genes:
        if g in _TIER2_ALIASES:
            gene_terms.append(_TIER2_ALIASES[g].pattern)
    gene_any_pat = re.compile(r"\b(?:" + "|".join(gene_terms) + r")\b", re.IGNORECASE)

    for sent in _split_sentences(text):
        if not any(p.search(sent) for p in GENE_CONTEXT_PATTERNS):
            continue
        # Skip sentences where the gene change is caused by a treatment/drug
        if _treatment_sentence(sent):
            continue

        words = sent.split()
        n = len(words)
        if n == 0:
            continue

        # Build character-start position for each word
        word_starts = []
        cursor = 0
        for w in words:
            idx = sent.find(w, cursor)
            word_starts.append(idx)
            cursor = idx + len(w)

        for pat in GENE_CONTEXT_PATTERNS:
            for m in pat.finditer(sent):
                # Map match span → word index range
                reg_words = [
                    wi for wi, ws in enumerate(word_starts)
                    if ws <= m.start() < ws + len(words[wi])
                    or m.start() <= ws < m.end()
                ]
                if not reg_words:
                    continue
                reg_lo, reg_hi = min(reg_words), max(reg_words)

                # Direct adjacency (1 word before or after the regulation span)
                for j in (reg_lo - 1, reg_hi + 1):
                    if 0 <= j < n and gene_any_pat.search(words[j]):
                        return True

                # 2 words away with exactly one linking word in between
                # e.g. "BDNF was upregulated", "expression of BDNF"
                for gene_j, mid_j in ((reg_lo - 2, reg_lo - 1),
                                       (reg_hi + 2, reg_hi + 1)):
                    if 0 <= gene_j < n and 0 <= mid_j < n:
                        mid = words[mid_j].lower().rstrip(".,;:")
                        if mid in _LINKING_WORDS and gene_any_pat.search(words[gene_j]):
                            return True

                # 3 words away with two intermediary words, at least one linking
                # e.g. "SNCA and BDNF were both expressed", "expression in the BDNF"
                for gene_j, mid1_j, mid2_j in (
                    (reg_lo - 3, reg_lo - 2, reg_lo - 1),
                    (reg_hi + 3, reg_hi + 1, reg_hi + 2),
                ):
                    if all(0 <= j < n for j in (gene_j, mid1_j, mid2_j)):
                        m1 = words[mid1_j].lower().rstrip(".,;:")
                        m2 = words[mid2_j].lower().rstrip(".,;:")
                        if (m1 in _LINKING_WORDS or m2 in _LINKING_WORDS) \
                                and gene_any_pat.search(words[gene_j]):
                            return True
    return False
def _split_sentences(text): return re.split(r"(?<=[.!?])\s+|\n", text)

# Section-header keywords that mark where the actual abstract body starts
_ABSTRACT_SECTION_RE = re.compile(
    r"(?m)^(?:BACKGROUND|OBJECTIVE[S]?|PURPOSE|INTRODUCTION|AIM[S]?|"
    r"CONTEXT|SUMMARY|METHODS?|RESULTS?|CONCLUSIONS?|ABSTRACT)[:\s]",
    re.IGNORECASE,
)

def _strip_pubmed_header(text):
    """Return only the abstract body, stripping the PubMed citation header
    (journal line, title, author names, affiliations).

    Strategy:
      1. If a section-header keyword (BACKGROUND:, OBJECTIVE:, etc.) is found,
         start from there.
      2. Otherwise find the end of the 'Author information:' block and take
         everything after the first blank line that follows.
      3. Fallback: return the full text unchanged.
    """
    m = _ABSTRACT_SECTION_RE.search(text)
    if m:
        return text[m.start():]
    # Find end of author information block
    ai = re.search(r"(?i)\bauthor information\b.*?\n(\n+)", text, re.DOTALL)
    if ai:
        return text[ai.end():]
    return text

def detect_mentioned_genes(text):
    """Detect any gene symbol mentioned in a biological context.

    Two paths:
    1. Alias-based: common names (e.g. "occludin", "tau") → canonical gene symbol
       via _TIER2_ALIASES, regardless of sentence context (the alias itself is specific).
    2. Symbol-based: regex finds any ALL-CAPS token; kept only when it appears in
       a sentence that also has biological vocabulary, and only if not in the
       non-gene abbreviation blacklist.
    """
    sentences = _split_sentences(text)
    found = set()

    # Path 1: alias / common-name detection
    for gene, alias_pat in _TIER2_ALIASES.items():
        if alias_pat.search(text):
            found.add(gene)

    # Path 2: symbol-level detection
    # Match on the ORIGINAL text — gene symbols are already ALL CAPS in papers;
    # lowercased anatomical terms (hippocampal, colonic) won't match.
    # Validate against HGNC approved symbols; resolve aliases to approved symbol.
    for sent in sentences:
        if not _BIO_CONTEXT_RE.search(sent):
            continue
        for m in _GENE_SYMBOL_RE.finditer(sent):
            token = m.group(1)
            if token in _NON_GENE_SYMBOLS:
                continue
            # Resolve: approved symbol or alias → approved symbol
            if token in HGNC_SYMBOLS:
                sym = token          # already the approved symbol (uppercase)
            elif token in _HGNC_ALIAS_TO_SYM:
                sym = _HGNC_ALIAS_TO_SYM[token].upper()
            else:
                continue
            if sym in BANNED_ACRONYM_CONTEXT:
                if any(re.search(banned, sent.lower()) for banned in BANNED_ACRONYM_CONTEXT[sym]):
                    continue
            found.add(sym)

    return found

# ─────────────────────────────────────────────────────────────────────────────
# CHRONOLOGICAL FETCH MODULE
# ─────────────────────────────────────────────────────────────────────────────
_API_KEY = os.environ.get("NCBI_API_KEY", None)
_SLEEP = 0.11 if _API_KEY else 0.34
if _API_KEY: Entrez.api_key = _API_KEY

DEBUG_YEAR_START = int(os.environ.get("DEBUG_YEAR_START", "2000"))   # full run = 2000
DEBUG_PER_YEAR   = os.environ.get("DEBUG_PER_YEAR")                  # None = no cap per year
if DEBUG_PER_YEAR is not None:
    DEBUG_PER_YEAR = int(DEBUG_PER_YEAR)

def fetch_abstracts(query, batch_size=500):
    current_year = datetime.datetime.now().year
    abstracts = []
    clean_query = re.sub(r"AND\s+\d{4}:\d{4}\[dp\]", "", query).strip()

    years = list(range(DEBUG_YEAR_START, current_year + 1))

    print(f"  Chronological splitting engine active ({DEBUG_YEAR_START} to {current_year}).")
    for year in years:
        yearly_query = f"({clean_query}) AND {year}[dp]"
        try:
            handle = Entrez.esearch(db="pubmed", term=yearly_query, retmax=100000)
            record = Entrez.read(handle)
            handle.close()
        except Exception as e:
            print(f"  Failed index initialization for year {year}: {e}")
            continue

        id_list = record["IdList"]
        yearly_total = len(id_list)
        if yearly_total == 0: continue

        if DEBUG_PER_YEAR is not None:
            id_list = id_list[:DEBUG_PER_YEAR]
            yearly_total = len(id_list)

        for start in range(0, yearly_total, batch_size):
            id_batch = id_list[start:start + batch_size]
            id_string = ",".join(id_batch)
            for attempt in range(3):
                try:
                    handle = Entrez.efetch(db="pubmed", id=id_string, rettype="abstract", retmode="text")
                    data = handle.read()
                    handle.close()
                    records = re.split(r"\n\n(?=\d+\.\s)", data)
                    abstracts.extend([r for r in records if len(r.strip()) > 50])
                    print(f"    Fetched {min(start + batch_size, yearly_total):>4}/{yearly_total} abstracts for year {year}", end="\r")
                    time.sleep(_SLEEP)
                    break
                except Exception as e:
                    time.sleep(2 ** (attempt + 1))
    print(f"\n  Done — Total of {len(abstracts)} abstracts aggregated.")
    return abstracts

# ── Multi-disorder detection: find which of our 10 disorders are mentioned ────
# Used to attribute a passing abstract to ALL relevant disorders, not just the
# one it was fetched under.
_DISORDER_KW_RE = {
    "Schizophrenia":    re.compile(r"\bschizophreni\w*\b|\bschizoaffective\b|\bSCZ\b", re.IGNORECASE),
    "Bipolar Disorder": re.compile(r"\bbipolar\b|\bmanic[\s\-]depress\w*\b|\bmanic\s+episode\b|\bBPD\b|\bBD\b", re.IGNORECASE),
    "Autism":           re.compile(r"\bautis\w*\b|\bASD\b|\bAsperger\b", re.IGNORECASE),
    "ADHD":             re.compile(r"\bADHD\b|\battention[\s\-]deficit\b", re.IGNORECASE),
    "Major Depression": re.compile(r"\bmajor\s+depress\w*\b|\bMDD\b|\bunipolar\s+depress\w*\b", re.IGNORECASE),
    "Anxiety":          re.compile(r"\banxiety\s+disorder\b|\bgeneralized\s+anxiety\b|\bsocial\s+anxiety\b|\bpanic\s+disorder\b", re.IGNORECASE),
    "PTSD":             re.compile(r"\bPTSD\b|\bpost[\s\-]?traumatic\s+stress\b|\bposttraumatic\s+stress\b", re.IGNORECASE),
    "OCD":              re.compile(r"\bOCD\b|\bobsessive[\s\-]compulsive\b", re.IGNORECASE),
    "Alzheimer's":      re.compile(r"\bAlzheimer\w*\b|\bAD\b(?!\s+libitum)", re.IGNORECASE),
    "Parkinson's":      re.compile(r"\bParkinson\w*\b|\bparkinsonism\b", re.IGNORECASE),
}

def detect_disorders_in_abstract(text):
    """Return a set of all disorders from DISORDERS that are mentioned in the abstract text."""
    return {d for d, pat in _DISORDER_KW_RE.items() if pat.search(text)}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN EXECUTION ENGINE
# ─────────────────────────────────────────────────────────────────────────────
results           = defaultdict(lambda: defaultdict(int))
tissue_results    = defaultdict(lambda: defaultdict(int))
gene_expr_counts  = defaultdict(int)
microbiome_counts = defaultdict(int)
gene_mention      = defaultdict(lambda: defaultdict(int))
disorder_genes    = defaultdict(set)
brain_gene_pool   = defaultdict(set)   # genes from brain-tissue + gene-context abstracts
gut_gene_pool     = defaultdict(set)   # genes from gut-tissue  + gene-context abstracts
rows              = []

_debug_dis = os.environ.get("DEBUG_DISORDERS")  # comma-sep disorder names -> restrict (testing)
_iter_disorders = ({d: DISORDERS[d] for d in [s.strip() for s in _debug_dis.split(",")] if d in DISORDERS}
                   if _debug_dis else DISORDERS)
for disorder, query in _iter_disorders.items():
    print(f"\nProcessing {disorder}…")
    fq = (f"{query} AND 2000:3000[dp] "
          f"NOT (clinical trial[pt] OR randomized controlled trial[pt] "
          f"OR meta-analysis[pt] OR review[pt])")
    abstracts = fetch_abstracts(fq)

    # ── Pass 1: build tissue + gene-context pools ─────────────────────────────
    # Three-phase split so the expensive LLM qualification call can be batched per
    # disorder (continuous batching) instead of firing sequentially per abstract.
    disorder_data = []
    seen_pmids = set()
    collected = []   # (item_id, abstract, body, tissue, has_gene_omics, has_micro, mentioned, prefilter_pass)

    # Phase A — cheap per-abstract regex features (no LLM). Tallies that are
    # independent of has_gene_context are updated here, exactly as before.
    for abstract in abstracts:
        # Deduplicate by PMID within this disorder
        pmid_m = re.search(r'\bPMID:\s*(\d+)', abstract)
        if pmid_m:
            pmid = pmid_m.group(1)
            if pmid in seen_pmids:
                continue
            seen_pmids.add(pmid)
            item_id = pmid
        else:
            # Stable id for records without a PMID line (cache key still valid).
            item_id = "b:" + hashlib.sha256(abstract.encode()).hexdigest()[:16]

        body = _strip_pubmed_header(abstract)

        # Skip long-term depression (LTD) papers misattributed to Major Depression.
        # (Disorder-specific cheap pre-filter; kept before the disorder-agnostic gate.)
        if _abstract_is_ltd_not_mdd(body, disorder):
            continue

        tissue          = detect_tissue(body)
        has_gene_omics  = detect_gene_expr(body)
        has_micro       = detect_microbiome(body)
        mentioned       = detect_mentioned_genes(body)
        prefilter_pass  = _regex_prefilter(body)   # cheap gate; LLM only sees passes (Option B)

        disorder_genes[disorder].update(mentioned)
        if has_gene_omics:  gene_expr_counts[disorder]  += 1
        if has_micro:       microbiome_counts[disorder] += 1
        for g in mentioned: gene_mention[disorder][g] += 1

        collected.append((item_id, abstract, body, tissue, has_gene_omics, has_micro, mentioned, prefilter_pass))

    # Phase B — batched qualification. LLM if enabled; otherwise (or on failure)
    # judgments stays None and Phase C falls back to the regex gate per abstract.
    # Option B: the LLM only judges pre-filter-passing abstracts (it replaces the
    # over-strict proximity heuristic on a bounded set), so the full fetch is not
    # LLM-judged — keeps the run to hours, not days.
    judgments = None
    if USE_LLM_QUALIFIER and _llm_qualify_fn is not None:
        items = [(iid, body, sorted(mentioned))
                 for (iid, _abs, body, _tis, _go, _mi, mentioned, pf) in collected if pf]
        print(f"[{disorder}] pre-filter {len(items)}/{len(collected)} -> LLM "
              f"(~{len(items)/19/60:.1f}h est @19 abstracts/min)")
        try:
            judgments = _llm_qualify_fn(items)
        except Exception as _e:
            print(f"[warn] LLM batch failed for {disorder} ({_e}); using regex gate.")
            judgments = None

    # Phase C — per-abstract gene-pool / tissue / regulation using the qualification
    # verdict. has_gene_context comes from the LLM (or the regex fallback); everything
    # downstream (tissue tally, per-gene tissue, the tissue-gene downgrade, multi-disorder
    # attribution, the appended tuple) runs identically to the original.
    for (item_id, abstract, body, tissue, has_gene_omics, has_micro, mentioned, prefilter_pass) in collected:
        if prefilter_pass and judgments is not None and item_id in judgments:
            has_gene_context = judgments[item_id].qualifies
        elif prefilter_pass:
            # LLM unavailable/failed for this disorder: fall back to the full
            # regex gate (incl. proximity) — identical to the pre-LLM behaviour.
            has_gene_context = detect_gene_context(body, mentioned)
        else:
            # Failed the cheap pre-filter (excluded confound / no disease model /
            # no expression signal) -> never qualifies; skips the LLM entirely.
            has_gene_context = False

        # Tissue tracking: split each tissue type into gene-context vs other (PRE-downgrade)
        if tissue != "Other":
            gc_tag = " (gene ctx)" if has_gene_context else " (other)"
            tissue_results[disorder][tissue + gc_tag] += 1
        else:
            tissue_results[disorder]["Other"] += 1

        # Per-gene tissue assignment (sentence-level for both brain and gut):
        # A gene goes to the brain pool only if it is mentioned near brain keywords.
        # A gene goes to the gut pool only if it is mentioned near gut keywords.
        genes_in_brain = set()
        genes_in_gut   = set()
        gene_regulation = {}   # gene → regulation direction for this abstract
        if has_gene_context:
            for gene in mentioned:
                gene_tissue = detect_tissue_for_gene(gene, body)
                if gene_tissue in ("Brain", "Brain & Gut"):
                    genes_in_brain.add(gene)
                    brain_gene_pool[disorder].add(gene)
                if gene_tissue in ("Gut", "Brain & Gut"):
                    genes_in_gut.add(gene)
                    gut_gene_pool[disorder].add(gene)
                gene_regulation[gene] = detect_gene_regulation(gene, body)

        is_brain_gc = bool(genes_in_brain)
        is_gut_gc   = bool(genes_in_gut)

        # Require at least one gene attributable to brain or gut tissue.
        # (Safety net: the LLM may qualify a serum-BDNF paper; the deterministic
        # per-gene tissue check finds no brain/gut gene and flips it False, so
        # fluid-only studies stay excluded exactly as under the regex gate.)
        if has_gene_context and not is_brain_gc and not is_gut_gc:
            has_gene_context = False

        # Multi-disorder attribution: if this abstract passes AND mentions other disorders,
        # attribute its genes to those disorders' pools too.
        if has_gene_context:
            other_disorders = detect_disorders_in_abstract(body) - {disorder}
            for other_d in other_disorders:
                if other_d in DISORDERS:
                    disorder_genes[other_d].update(mentioned)
                    for g in mentioned:
                        gene_mention[other_d][g] += 1
                    for g in genes_in_brain:
                        brain_gene_pool[other_d].add(g)
                    for g in genes_in_gut:
                        gut_gene_pool[other_d].add(g)

        disorder_data.append(
            (abstract, body, tissue, has_gene_context, has_micro,
             mentioned, is_brain_gc, is_gut_gc, genes_in_brain, genes_in_gut,
             gene_regulation)
        )

    # ── Shared pool: genes in BOTH brain and gut gene-context abstracts ───────
    shared_pool = brain_gene_pool[disorder] & gut_gene_pool[disorder]

    # ── Pass 2: annotate rows using the completed pools ───────────────────────
    for abstract, body, tissue, has_gc, has_micro, mentioned, is_brain_gc, is_gut_gc, genes_in_brain, genes_in_gut, gene_regulation \
            in disorder_data:
        # Use per-gene tissue sets directly for brain/gut found
        brain_found  = sorted(genes_in_brain) if is_brain_gc else []
        gut_found    = sorted(genes_in_gut)   if is_gut_gc   else []
        # Shared: genes that appear in brain context in this abstract AND gut context elsewhere (or vice versa)
        shared_found = sorted((genes_in_brain | genes_in_gut) & shared_pool) if (is_brain_gc or is_gut_gc) else []

        # Build base row dict
        base_row = {
            "abstract":           body,
            "tissue":             tissue,
            "gene_context":       has_gc,
            "microbiome":         has_micro,
            "genes_found":        ", ".join(sorted(mentioned)),
            "brain_genes_found":  ", ".join(brain_found),
            "gut_genes_found":    ", ".join(gut_found),
            "shared_genes_found": ", ".join(shared_found),
            "gene_regulation":    str({g: gene_regulation.get(g, "Unknown") for g in mentioned}) if gene_regulation else "",
        }
        rows.append({"disorder": disorder, **base_row})

        # Multi-disorder attribution: add duplicate rows for other detected disorders
        if has_gc:
            other_disorders = detect_disorders_in_abstract(body) - {disorder}
            for other_d in other_disorders:
                if other_d in DISORDERS:
                    rows.append({"disorder": other_d, **base_row})

df = pd.DataFrame(rows)
df.to_csv("model_organisms_dataset.csv", index=False)

# ── Save three text files — one per gene category ─────────────────────────────
SEP = "=" * 60

def save_as_text(df_sub, gene_col, out_path):
    lines = []
    for idx, (_, row) in enumerate(df_sub.iterrows()):
        genes    = row[gene_col]   if pd.notna(row.get(gene_col, "")) else ""
        abstract = row["abstract"] if pd.notna(row["abstract"])       else ""
        lines.append(f"[{idx}, {row['disorder']}, {genes}, [{abstract}]]")
        lines.append(SEP)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Saved {out_path}  ({len(df_sub):,} abstracts)")

df_brain = df[df["brain_genes_found"].str.len() > 0].copy()
save_as_text(df_brain, "brain_genes_found", "abstracts_brain_genes_context.txt")

df_gut = df[df["gut_genes_found"].str.len() > 0].copy()
save_as_text(df_gut, "gut_genes_found", "abstracts_gut_genes_context.txt")

df_both = df[df["shared_genes_found"].str.len() > 0].copy()
save_as_text(df_both, "shared_genes_found", "abstracts_shared_genes_context.txt")

# ── All passing abstracts (gene_context=True) in one CSV ──────────────────────
# Extract publication year from abstract copyright line
def _extract_year(text):
    m = re.search(r'(?:©|\(C\)|Copyright)\s+(?:\d{4}\s+)?(?:[A-Za-z\s,\.]+\s+)?(\d{4})', text)
    if m:
        return int(m.group(1))
    m2 = re.search(r'\b(20[0-2]\d|199\d)\b', text[-200:])
    if m2:
        return int(m2.group(1))
    return None

df_passing = df[df["gene_context"] == True].copy()
df_passing["year"] = df_passing["abstract"].apply(_extract_year)
pass_cols = ["disorder", "year", "tissue", "brain_genes_found", "gut_genes_found",
             "shared_genes_found", "gene_regulation", "abstract"]
df_passing[pass_cols].to_csv("passing_abstracts.csv", index=False)
print(f"Saved passing_abstracts.csv  ({len(df_passing):,} abstracts that passed all filters)")

# ── Summary file: abstract counts per disorder in the shared file ──────────────
summary_lines = ["Shared Gene-Context Abstracts — Count per Disorder", "=" * 50]
counts = df_both["disorder"].value_counts().reindex(list(DISORDERS.keys()), fill_value=0)
for disorder, count in counts.items():
    summary_lines.append(f"{disorder:<20} {count:>6} abstracts")
summary_lines.append("=" * 50)
summary_lines.append(f"{'TOTAL':<20} {counts.sum():>6} abstracts")
with open("shared_genes_summary.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(summary_lines))
print(f"Saved shared_genes_summary.txt")

# ── Per-disorder gene sets (data-driven from tissue + gene-context abstracts) ─
disorder_brain_only_genes = {d: brain_gene_pool[d] - gut_gene_pool[d] for d in DISORDERS}
disorder_gut_only_genes   = {d: gut_gene_pool[d]   - brain_gene_pool[d] for d in DISORDERS}
disorder_shared_genes     = {d: brain_gene_pool[d] & gut_gene_pool[d]   for d in DISORDERS}




# ─────────────────────────────────────────────────────────────────────────────
# TO-DO DISORDER PAIRING & SHARED GENE REPORTS
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*50 + "\nROADMAP METRICS AND PAIR SELECTIONS\n" + "="*50)

all_disorders = list(DISORDERS.keys())

def _safe_intersect(sets):
    filled = [s for s in sets if s]
    return set.intersection(*filled) if filled else set()

shared_brain_only = _safe_intersect([disorder_brain_only_genes[d] for d in all_disorders])
shared_gut_only   = _safe_intersect([disorder_gut_only_genes[d]   for d in all_disorders])
shared_both       = _safe_intersect([disorder_shared_genes[d]     for d in all_disorders])

print(f"Shared Brain-only Genes across ALL disorders: {sorted(shared_brain_only) or '(none)'}")
print(f"Shared Gut-only Genes across ALL disorders:   {sorted(shared_gut_only)   or '(none)'}")
print(f"Shared Brain+Gut Genes across ALL disorders:  {sorted(shared_both)       or '(none)'}")

pairing_scores = []
processed_pairs = set()

for d1 in DISORDERS:
    for d2 in DISORDERS:
        if d1 != d2 and (d2, d1) not in processed_pairs:
            shared_set = disorder_genes[d1].intersection(disorder_genes[d2])
            pairing_scores.append((d1, d2, len(shared_set), sorted(list(shared_set))))
            processed_pairs.add((d1, d2))

pairing_scores.sort(key=lambda x: -x[2])

print("\n The pairs of mental disorders sharing the most genes:")
for i, (dis1, dis2, score, shared_list) in enumerate(pairing_scores):
    print(f"  {i+1}. {dis1} + {dis2} -> Shared Count: {score} genes")
    print(f"     Genes: {', '.join(shared_list[:12])}...")


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING HELPER
# ─────────────────────────────────────────────────────────────────────────────
def venn_label(gene_list, count_dict, max_show=7):
    shown = gene_list[:max_show]
    lines = [f"{g} ({count_dict.get(g, 0)})" for g in shown]
    extra = len(gene_list) - max_show
    return "\n".join(lines) + (f"\n…+{extra} more" if extra > 0 else "")


# ─────────────────────────────────────────────────────────────────────────────
# PER-DISORDER PLOTTING GENERATOR (Species pie + Tissue pie + Venn)
# ─────────────────────────────────────────────────────────────────────────────
rodent_names = set(RODENT_SPECIES.keys()) | {"Other Rodents"}
blue_shades  = plt.cm.Blues([0.9, 0.75, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1])
other_colors = list(plt.cm.Set2.colors)

print("\nGenerating Diagnostic Plots...")
for disorder in DISORDERS:
    data = results[disorder]
    if not data:
        continue

    fig, (ax_ti, ax_venn) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(disorder, fontsize=16, fontweight="bold", y=1.01)

    # Panel 1 — Tissue pie plot (gene-context breakdown per tissue type)
    td = tissue_results[disorder]
    t_order = [
        "Brain (gene ctx)", "Brain (other)",
        "Gut (gene ctx)",   "Gut (other)",
        "Brain & Gut (gene ctx)", "Brain & Gut (other)",
    ]
    t_cmap = {
        "Brain (gene ctx)":       "#1f4e79",
        "Brain (other)":          "#9dc3e6",
        "Gut (gene ctx)":         "#833c00",
        "Gut (other)":            "#f4b183",
        "Brain & Gut (gene ctx)": "#375623",
        "Brain & Gut (other)":    "#a9d18e",
    }
    tl, ts, tc = [], [], []
    for t in t_order:
        if td.get(t, 0) > 0:
            tl.append(t); ts.append(td[t]); tc.append(t_cmap[t])
    if ts:
        w2, _ = ax_ti.pie(ts, colors=tc, startangle=140)
        ax_ti.legend(w2, [f"{l} (n={s})" for l, s in zip(tl, ts)],
                     title="Research Focus\n(dark = gene ctx)", loc="center left",
                     bbox_to_anchor=(1, 0.5), fontsize=8)
    ax_ti.set_title(f"Brain vs Gut Research Focus\n{disorder} (2000–present)")
    ax_ti.set_aspect("equal")

    # Panel 3 — Venn
    # Left   : genes found only in brain gene-context abstracts
    # Right  : genes found only in gut gene-context abstracts
    # Overlap: genes found in both brain AND gut gene-context abstracts
    dc = gene_mention[disorder]
    bp = brain_gene_pool[disorder]
    gp = gut_gene_pool[disorder]
    sp = bp & gp
    fb = sorted(bp - sp, key=lambda g: -dc.get(g, 0))
    fg = sorted(gp - sp, key=lambda g: -dc.get(g, 0))
    fo = sorted(sp,       key=lambda g: -dc.get(g, 0))

    vd = venn2(
        subsets=(len(fb) or 1, len(fg) or 1, len(fo) or 1),
        set_labels=("Brain-only\nGenes", "Gut-only\nGenes"),
        ax=ax_venn,
        set_colors=("#4C72B0", "#DD8452"),
        alpha=0.55,
    )
    for rid, gl in [("10", fb), ("01", fg), ("11", fo)]:
        lbl = vd.get_label_by_id(rid)
        if lbl:
            shown = gl[:7]
            text  = "\n".join(f"{g}({dc.get(g,0)})" for g in shown)
            if len(gl) > 7: text += f"\n…+{len(gl)-7}"
            lbl.set_text(text)
            lbl.set_fontsize(6)
    ax_venn.set_title(
        f"Genes: Brain-only / Gut-only / Shared\n{disorder}\n"
        f"(from gene-context abstracts, by tissue)",
        fontsize=9)

    plt.tight_layout()
    fname = ("model_organisms_"
             + disorder.replace(" ", "_").replace("/", "_").replace("'", "") + ".png")
    plt.savefig(fname, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"Saved visualization file: {fname}")

    # ── Per-disorder gene lists (console) ────────────────────────────────────
    # Brain list = brain-only genes + shared genes (both belong to brain)
    # Gut list   = gut-only genes   + shared genes (both belong to gut)
    # Shared     = genes in both lists (shown separately for clarity)
    fb_full = sorted(fb + fo, key=lambda g: -dc.get(g, 0))
    fg_full = sorted(fg + fo, key=lambda g: -dc.get(g, 0))

    print(f"\n   {disorder}: Brain genes (brain-only + shared)")
    for g in fb_full or ["(none detected)"]:
        tag = " [shared]" if g in set(fo) else ""
        print(f"    {g:<12}  {dc.get(g,0):>5} mentions{tag}")

    print(f"   {disorder}: Gut genes (gut-only + shared)")
    for g in fg_full or ["(none detected)"]:
        tag = " [shared]" if g in set(fo) else ""
        print(f"    {g:<12}  {dc.get(g,0):>5} mentions{tag}")

    print(f"   {disorder}: Brain & Gut shared genes")
    for g in fo or ["(none detected)"]:
        print(f"    {g:<12}  {dc.get(g,0):>5} mentions")