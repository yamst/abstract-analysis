import re, time, datetime, os, ssl, socket, hashlib
import urllib.request
from Bio import Entrez
from collections import defaultdict

socket.setdefaulttimeout(60)
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "dasfour@outlook.com"

DEBUG_DISORDER = "Schizophrenia"
DEBUG_QUERY    = "schizophrenia AND gene expression"
DEBUG_N        = 20          # number of abstracts to fetch
DEBUG_YEAR     = 2018        # single year to keep the fetch tiny

BRAIN_GENES = list(dict.fromkeys([
    "RBFOX3","MAP2","TUBB3","SYP","SYT1","DLG4","GFAP","MBP","AIF1","ENO2",
    "NEFH","NEFL","SNAP25","TH","GAD1","SLC17A7","PVALB","SST","OLIG2","S100B",
    "DCX","STMN2","ASCL1","PROX1","NRGN","CAMK2A","SLC17A6","GAD2","SLC32A1",
    "FOXP2","ALDH1L1","SLC1A2","AQP4","MOG","PLP1","CNP","TMEM119","P2RY12","CX3CR1",
    "BDNF","NGF","NTF3","NTRK2","GDNF","NPY","VIP","GAL","ADCYAP1","CALCA",
    "TAC1","PENK","SST","NOS1","CRH","NR3C1","TPH2","SLC6A4","HTR1A","HTR2B",
    "HTR3A","HTR4","DRD2","TLR4","TNF","IL6","IL1B","HMGB1","CHAT","IDO1",
    "CLOCK","PER2","ARNTL","NR1D1","LEPR","CCK","AGRP","SYN1","NCAM1",
    "SNCA","LRRK2","PRKN","PINK1","BECN1","SQSTM1","ACE2","SHANK3","CHD8",
    "CNTNAP2","NRXN1","NLGN3","SYNGAP1","DYRK1A","SCN2A","FOXP1","TCF4","ANK2",
    "RET","PHOX2B","TJP1","FFAR2","VIPR1","GLP1R","GLP2R","ARNTL",
]))
GUT_GENES = list(dict.fromkeys([
    "VIL1","CDX2","LGR5","MUC2","ATOH1","DEFA5","CHGA","GIP","SLC5A1","FABP2",
    "SI","EPCAM","CLDN3","OCLN","OLFM4","REG3A","AQP3","SLC26A3","GLP1R","PYY",
    "ASCL2","SMOC2","DEFA6","LYZ","CHGB","NEUROG3","PAX4","INSM1","GCG","SCT",
    "GHRL","TPH1","CA2","AQP8","SLC26A2","CLDN4","CLDN7","CDH1","MUC5B","TFF3",
    "FCGBP","ANPEP","SLC15A1","APOA4","APOB","ALPI","SLC2A5","LEP",
    "BDNF","NGF","NTF3","NTRK2","GDNF","NPY","VIP","GAL","ADCYAP1","CALCA",
    "TAC1","PENK","SST","NOS1","CRH","NR3C1","TPH1","SLC6A4","HTR1A","HTR2B",
    "HTR3A","HTR4","DRD2","TLR4","TNF","IL6","IL1B","HMGB1","CHAT","IDO1",
    "CLOCK","PER2","ARNTL","NR1D1","LEPR","CCK","AGRP","SYN1","NCAM1",
    "SNCA","LRRK2","PRKN","PINK1","BECN1","SQSTM1","ACE2","SHANK3","CHD8",
    "CNTNAP2","NRXN1","NLGN3","SYNGAP1","DYRK1A","SCN2A","FOXP1","TCF4","ANK2",
    "RET","PHOX2B","TJP1","FFAR2","VIPR1","GLP1R","GLP2R",
]))
ALL_CANDIDATE_GENES = list(dict.fromkeys(BRAIN_GENES + GUT_GENES))

TREATMENT_CAUSAL_RE_OLD = re.compile(
    r"\btreated\s+with\b"
    r"|\btreatment\s+with\b"
    r"|\badministration\s+of\b"
    r"|\b(reversed|restored|normalized|attenuated|rescued|ameliorated)\s+by\b"
    r"|\bin\s+response\s+to\b"
    r"|\bfollowing\s+\w+\s*(treatment|administration|injection|infusion)\b"
    r"|\bafter\s+\w+\s*(treatment|administration|injection)\b"
    r"|\binduced\s+by\b"
    r"|\b\d+[\s]*(mg/kg|μg|μM|nM|mg/ml|nmol)\b",
    re.IGNORECASE
)

_DRUG_KW = (
    r"drug[s]?|medication[s]?|antidepressant[s]?|antipsychotic[s]?|"
    r"anxiolytic[s]?|stimulant[s]?|vehicle|placebo|"
    r"haloperidol|clozapine|risperidone|olanzapine|quetiapine|"
    r"aripiprazole|ziprasidone|amisulpride|paliperidone|lurasidone|"
    r"chlorpromazine|fluphenazine|perphenazine|thioridazine|"
    r"fluoxetine|sertraline|paroxetine|escitalopram|citalopram|"
    r"venlafaxine|duloxetine|mirtazapine|bupropion|trazodone|"
    r"imipramine|amitriptyline|nortriptyline|clomipramine|"
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
    r"bisphenol[\s\-]?a|\bbpa\b|bisphenol[\s\-]?s|\bbps\b|"
    r"neuroleptic[s]?|"
    r"n[\s\-]?acetylcysteine|\bnac\b|"
    r"polyphenol[s]?|theogallin|resveratrol|curcumin|quercetin|"
    r"sigma[\s\-]?\d+\s+receptor\s+(agonist|antagonist)|"
    r"\bPRE084\b|\bDTG\b|\bBD1047\b|"
    r"herbal\s+(medicine|extract|formula|remedy|preparation)|"
    r"traditional\s+chinese\s+medicine|\bTCM\b|"
    r"oral\s+liquid|oral\s+decoction|"
    r"plant\s+extract|phytochem\w+|"
    r"aloe[\s\-]vera|"
    r"\bintragastric\s+administration\b|\bgavage\s+(?:administration|group)\b|"
    r"\btreatment\s+group\s+(?:received|was\s+given|were\s+given|was\s+administered)\b|"
    r"\bPills?\s+(?:were|was)\s+given\b|"
    r"ascorbic\s+acid|"
    r"propionic\s+acid|\bshort[\s\-]chain\s+fatty\s+acid[s]?\b|"
    r"wheat\s+(malt|germ)\s+extract|food[\s\-]derived\s+(source|extract)|"
    r"nutritional\s+strateg\w+|dietary\s+supplement\w*|"
    r"memantine|donepezil|galantamine|rivastigmine|aducanumab|lecanemab|"
    r"levodopa|l[\s\-]?dopa|carbidopa|pramipexole|ropinirole|selegiline|rasagiline|"
    r"retinoic\s+acid|\bngf\s+treatment\b|"
    r"d[\s\-]?galactose|hydrogen\s+peroxide|\bH2O2\b|"
    r"ceftriaxone|nanoparticl\w+|"
    r"vitamin\s+[ABCDE]\s+(supplement|supplementation|treatment|deficiency\s+rescue)|"
    r"lps|lipopolysaccharide(?!\s+receptor)"
)

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
    r"trimethyltin|tmt|"
    r"colchicine|"
    r"paraquat|maneb"
)

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
_MI_HYPHEN_RE = re.compile(
    r"\b(?:" + _MODEL_INDUCER_KW + r")[\s\-](?:treated|induced|challenged)\b",
    re.IGNORECASE,
)
_MI_BARE_RE = re.compile(
    r"\b(?:" + _MODEL_INDUCER_KW + r")\b.{0,200}"
    r"(?:stimulat\w+|induc\w+|trigger\w+|promot\w+|caus\w+)\s+\w[\w\s]{0,30}"
    r"(?:gene|mRNA|expression|transcript)",
    re.IGNORECASE | re.DOTALL,
)

def _abstract_uses_model_inducer(text):
    return bool(
        _MI_ADMIN_RE.search(text)
        or _MI_ADMIN_RE2.search(text)
        or _MI_HYPHEN_RE.search(text)
        or _MI_BARE_RE.search(text)
    )

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
_TD_HYPHEN_RE = re.compile(
    r"\b(?:" + _DRUG_KW + r")[\s\-](?:treated|induced|mediated)\b",
    re.IGNORECASE,
)

def _abstract_uses_therapeutic_drug(text):
    if _TD_ADMIN_RE.search(text) or _TD_ADMIN_RE2.search(text):
        return True
    if _DOSAGE_RE.search(text):
        return True
    if _TD_HYPHEN_RE.search(text):
        return True
    if _TD_DRUG_ANYWHERE_RE.search(text) and _TD_ADMIN_ANYWHERE_RE.search(text):
        return True
    return False

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
    r"\bknockdown\s+(of|in|using|via|with|by)\b|"
    r"\boverexpression\s+vector\b|"
    # iPSC removed here — patient-derived iPSC is valid, caught by positive filter
    # r"\binduced\s+pluripotent\s+stem\s+cell[s]?\b|\biPSC[s]?\b|"
    r"\b(?:reprogramming|reprogrammed)\s+(?:somatic\s+)?cells?\b|"
    r"\b(?:cord\s+blood|urine|saliva|hair\s+follicle)\s+(?:cells?\s+)?(?:for|as)\s+(?:a\s+)?(?:source|reprogramming)\b|"
    r"\biPSC\s+(?:derivation|generation|reprogramming|technology)\b|"
    r"\bgeneration\s+of\s+iPSC\b|"
    r"\btreated\s+(?:the\s+)?(?:cells?|neurons?|cultures?|slices?|astrocytes?)\s+with\b|"
    r"\bstereotaxic\s+(injection|surgery|implant)\b|"
    r"\bintracranial\s+(injection|infusion|implant)\b|"
    r"\bretrograde\s+(tracing|virus|vector)\b|"
    r"\banterograde\s+(tracing|virus|vector)\b|"
    # ── Behavioural / psychological stress protocols ───────────────────────────
    r"\bforced\s+swim\b|\btail\s+suspension\b|"
    r"\b(im)?mobilization\s+stress\b|\brestraint\s+stress\b|"
    r"\bchronic\s+(mild|unpredictable|variable|social|restraint|defeat)\s+stress\b|"
    r"\bsocial\s+(defeat|isolation\s+stress)\b|\bpredator\s+(stress|exposure)\b|"
    r"\bearly\s+(life\s+)?stress\b|\bneonatal\s+(stress|separation|isolation)\b|"
    r"\bmaternal\s+(separation|deprivation)\b|"
    r"\bfear\s+conditioning\b|\bconditioned\s+fear\b|"
    r"\bassociative\s+learning\b|\bpavlovian\b|"
    # ── Dietary / metabolic interventions ─────────────────────────────────────
    r"\bhigh[\s\-]fat\s+diet\b|\bketogenic\s+diet\b|\bcaloric\s+restriction\b|"
    r"\bdietary\s+(intervention|restriction|supplementation|deprivation)\b|"
    r"\bdietary\s+(pufa|omega[\s\-]?3|fatty\s+acid)\s+deprivation\b|"
    r"\bn[\s\-]?3\s+pufa\s+deprivation\b|"
    r"\bretraction\s+in\b|"
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
    r"\btransgenic\s+(livestock|animal|overexpression|construct)\b|"
    r"\bsite[\s\-]specific\s+integration\b|"
    r"\btransgene\s+(integration|expression|insert)\b|"
    # ── Irradiation ───────────────────────────────────────────────────────────
    r"\birradiation\b|\bradiation\s+(exposure|treatment)\b",
    re.IGNORECASE,
)

def _abstract_uses_nondrug_intervention(text):
    return bool(_NONDRUG_INTERVENTION_RE.search(text))

_DISEASE_MODEL_RE = re.compile(
    r"\bpatients?\s+with\b|\bdiagnosed\s+with\b|"
    r"\bsubjects?\s+with\b|\bindividuals?\s+with\b|"
    r"\bpostmortem\b|\bpost[\s\-]mortem\b|\bautopsy\b|\bbrain\s+bank\b|"
    r"\bbrain\s+tissue\s+from\b|\btissue\s+samples?\s+from\b|"
    r"\bpatient[\s\-](?:derived|specific)\b|"
    r"\bfrom\s+(?:patients?|affected\s+individuals?|diagnosed\s+(?:individuals?|subjects?))\b|"
    r"\bcase[\s\-]control\b|\bclinical\s+sample[s]?\b|"
    r"\bschizophrenia\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bbipolar\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bautism\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bADHD\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bdepression\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\banxiety\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bPTSD\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bOCD\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bAlzheimer\w*\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\bParkinson\w*\s+(?:patients?|subjects?|group|brain|cohort)\b|"
    r"\b(?:mouse|rat|animal|primate|rodent|murine)\s+model\s+of\b|"
    r"\bmodel\s+(?:of|for)\s+\w[\w\s]{0,30}(?:disorder|disease|syndrome)\b|"
    r"\b(?:disorder|disease)[\s\-]model\b|"
    r"\b(?:knock[\s\-]?out|knock[\s\-]?in)\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\bKO\s+(?:mice|rats|mouse|rat)\b|"
    r"\btransgenic\s+(?:mice|rats|mouse|rat|animals?|model|line)\b|"
    r"\bmutant\s+(?:mice|rats|mouse|rat)\b|"
    r"\bhaploinsuffici\w+\s+(?:mice|rats|mouse|rat|animals?)\b|"
    r"\bAPP[\s\-/]?PS\d\b|\b5xFAD\b|\b3xTg[\s\-]?AD\b|"
    r"\bLRRK2\s+(?:G2019S|R1441|I2020T|Y1699C)\b|"
    r"\bShank\d[\s\-](?:KO|knockout|deficient|mutant)\b|"
    r"\b6[\s\-]?(?:ohda|hydroxydopamine)\b|\bmptp\b|"
    r"\bstreptozotocin\b|\bstz\s+(?:mice|rats|model)\b|"
    r"\bpoly[\s\-]?i[\s\-:]?c\b|\bmk[\s\-]?801\b|"
    r"\bscopolamine\b|\bcuprizone\b|\bkainic\s+acid\b|"
    r"\bquinolinic\s+acid\b|\brotenone\b|\bparaquat\b|\breserpine\b|"
    r"\bokadaic\s+acid\b|\blactacystin\b|"
    r"\biPSC[\s\-]?derived\b|\bpatient[\s\-]derived\s+(?:neurons?|cells?|iPSC)\b|"
    r"\bpatient[\s\-]specific\s+(?:neurons?|cells?|iPSC)\b|"
    r"\biPSC\s+from\s+(?:patients?|individuals?|subjects?)\b|"
    r"\bdisease[\s\-]specific\s+(?:neurons?|cells?|iPSC)\b",
    re.IGNORECASE,
)

def _has_disease_model(text):
    return bool(_DISEASE_MODEL_RE.search(text))

_BEHAVIORAL_AIM_RE = re.compile(
    r"\baim(?:ed|s)?\s+(?:to\s+)?(?:was\s+to\s+)?(?:investigate|examine|assess|evaluate|study|"
    r"explore|determine|test|characterize|measure|quantify)\s+.{0,80}"
    r"(?:behavioral?|locomotor|cognitive|social\s+(?:behavior|interaction|memory)|"
    r"memory|learning|anxiety[\s\-]like|depression[\s\-]like|fear|attention)\s+"
    r"(?:deficits?|impairment|performance|phenotype|function|outcome|changes?|capacit\w+|abilit\w+)\b|"
    r"\bwe\s+(?:investigated|examined|assessed|evaluated|studied|explored|measured|tested|"
    r"characterized)\s+.{0,60}"
    r"(?:behavioral?|locomotor|cognitive|social)\s+"
    r"(?:deficits?|impairment|performance|phenotype|function|outcome|changes?)\b|"
    r"\b(?:behavioral?|locomotor|cognitive|social)\s+"
    r"(?:deficits?|impairment|performance|changes?|phenotype)\s+"
    r"(?:were|was)\s+(?:the\s+)?(?:primary|main|principal|key|major)\s+"
    r"(?:endpoint|outcome|finding|result|readout)\b|"
    r"\bbehavioral?\s+(?:deficits?|impairment|changes?|outcomes?)\s+were\s+"
    r"(?:observed|found|detected|demonstrated|shown|noted|documented)\b",
    re.IGNORECASE | re.DOTALL,
)

_EXPRESSION_MEASUREMENT_RE = re.compile(
    r"\bmRNA\b|"
    r"\bRT[\s\-]?PCR\b|\bqPCR\b|q[\s\-]?RT[\s\-]?PCR\b|"
    r"\bmicroarray\b|\bRNA[\s\-]?seq\b|\bscRNA[\s\-]?seq\b|"
    r"\btranscriptom\w+\b|"
    r"\bin\s+situ\s+hybridization\b|\bISH\b|"
    r"\bWestern\s+blot\b|\bimmunoblot\b|"
    r"\bimmunohistochem\w+\b|\bimmunofluo\w+\b|"
    r"\bexpression\s+profil\w+\b|"
    r"\b(?:up|down)[\s\-]?regulat\w+\b|"
    r"\bdifferentially\s+expressed\b|\bdifferential\s+(?:gene\s+)?expression\b|"
    r"\boverexpressed\b|\bunderexpressed\b|"
    r"\bexpression\s+(?:was|were|is|are)\s+\w*\s*"
    r"(?:increased|decreased|reduced|elevated|altered|changed|upregulated|"
    r"downregulated|measured|quantified|higher|lower|suppressed|enhanced|"
    r"significantly\s+\w+)\b|"
    r"\b(?:increased|decreased|reduced|elevated|altered|changed|higher|lower|"
    r"enhanced|suppressed|significant)\s+(?:\w+\s+){0,3}expression\b|"
    r"\bexpression\s+levels?\s+(?:were|was|are|is)?\s*"
    r"(?:increased|decreased|reduced|elevated|altered|changed|higher|lower|"
    r"significantly\s+\w+)\b|"
    r"\bmRNA\s+levels?\s+(?:were|was|are|is)?\s*"
    r"(?:increased|decreased|reduced|elevated|altered|changed|higher|lower)\b|"
    r"\bgene\s+expression\s+(?:analysis|profil\w+|data|levels?|changes?|pattern|studi\w+)\b",
    re.IGNORECASE,
)

def _abstract_is_behavioral_primary(text):
    if not _BEHAVIORAL_AIM_RE.search(text):
        return False
    return not _EXPRESSION_MEASUREMENT_RE.search(text)

_GENETIC_ASSOC_RE = re.compile(
    r"\b(?:polymorphism[s]?|SNP[s]?|genotype[s]?|haplotype[s]?|allele[s]?|variant[s]?|"
    r"copy\s+number\s+variation[s]?|\bCNV[s]?\b)\s+.{0,60}"
    r"(?:(?:is|are|was|were)\s+(?:not\s+)?(?:significantly\s+)?associated\s+with|"
    r"(?:significantly\s+)?predict[s]?\s+risk|confer[s]?\s+(?:increased\s+)?risk)\b|"
    r"\bgenome[\s\-]wide\s+association\b|\bGWAS\b|"
    r"\blinkage\s+(?:disequilibrium|analysis|study)\b|\bLD\s+block\b|"
    r"\bassociation\s+(?:study|analysis|between\s+\w+\s+(?:polymorphism|SNP|genotype))\b|"
    r"\bno\s+(?:significant\s+)?association\s+(?:was\s+found\s+)?between\s+\w[\w\s]{0,30}"
    r"(?:polymorphism|SNP|genotype|variant|allele)\b|"
    r"\brisk\s+(?:allele[s]?|variant[s]?|genotype[s]?|haplotype[s]?)\b|"
    r"\bsusceptibility\s+(?:allele[s]?|variant[s]?|locus|loci)\b|"
    r"\b(?:disease|disorder)\s+risk\s+(?:allele|variant|SNP|genotype)\b",
    re.IGNORECASE | re.DOTALL,
)

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
    return False

_EPIGENETICS_RE = re.compile(
    r"\bepigeneti[ck]\w*\b|\bepigenome\b|\bepigenomic\w*\b|"
    r"\b(de|un|re|hyper|hypo)?methylat\w+\b|"
    r"\bDNA\s+methylation\b|\bCpG\b|"
    r"\b(de)?acetylat\w+\b|"
    r"\bubiquitinat\w+\b|\bubiquitylat\w+\b|"
    r"\bhistone\s+(modification|mark|acetyl|methyl|ubiquitin|variant)\w*\b|"
    r"\bchromatin\s+(remodel\w+|modifi\w+|accessib\w+|structure)\b|"
    r"\bchromatin\s+immunoprecipitation\b|\bChIP\b|"
    r"\bmiRNA[s\-]?\b|\bmicro[\s\-]?RNA[s]?\b|\bmiR[\s\-]\d+\b|"
    r"\bsRNA[s]?\b|\bsmall[\s\-]RNA[s]?\b|"
    r"\bnon[\s\-]coding[\s\-]RNA[s]?\b|\bncRNA[s]?\b|"
    r"\blncRNA[s]?\b|\blong[\s\-]non[\s\-]coding[\s\-]RNA[s]?\b|"
    r"\bcircRNA[s]?\b|\bpiRNA[s]?\b|"
    r"\bSUMOylation\b|\bsumoylat\w+\b",
    re.IGNORECASE,
)

def _abstract_uses_epigenetics(text):
    return bool(_EPIGENETICS_RE.search(text))

_DRUG_TREATMENT_NOUN_RE = re.compile(
    r"\b(?:" + _DRUG_KW + r")\s+treatment\b|"
    r"\btreatment\s+with\s+(?:" + _DRUG_KW + r")\b|"
    r"\bduring\s+(?:" + _DRUG_KW + r")\s+treatment\b",
    re.IGNORECASE,
)
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
    r"\btreatment\s+with\s+[\w\s,\-\+\(\)]{2,60}"
    r"(significantly\s+)?(increased?|decreased?|reversed?|improved?|restored?|attenuated?|"
    r"enhanced?|reduced?|elevated?|normalized?|rescued?|prevented?|altered?|modulated?)\b|"
    r"\b(treated|treatment)\s+with\s+\w[\w\s,\-]+(increased?|decreased?|altered?|"
    r"upregulated?|downregulated?|modulated?)\s+(gene\s+|mRNA\s+|protein\s+)?expression\b|"
    r"\bchallenged\s+with\s+\w[\w\s\-,]+(?:inhibitor|agonist|antagonist|compound|drug|agent)\b|"
    r"\btreatment\s+with\s+(the\s+)?recombinant\b|"
    r"\btreated\s+with\s+(the\s+)?recombinant\b|"
    r"\brecombinant\s+\w[\w\-]+\s+(protein|peptide|factor|cytokine)\b|"
    r"\b(?:" + _DRUG_KW + r")\s+(significantly\s+)?"
    r"(increased?|decreased?|upregulated?|downregulated?|altered?|modulated?|"
    r"induced?|suppressed?|enhanced?|reduced?)\s+(the\s+)?"
    r"(gene\s+|mRNA\s+|protein\s+)?expression\b",
    re.IGNORECASE | re.DOTALL,
)

def _abstract_attributes_expression_to_treatment(text):
    return bool(
        _TREATMENT_ATTRIBUTION_RE.search(text)
        or _THERAPEUTIC_EFFECTS_RE.search(text)
        or _DRUG_TREATMENT_NOUN_RE.search(text)
        or _LI_TREATMENT_RE.search(text)
    )

TREATMENT_CAUSAL_RE_NEW = re.compile(
    # Direct administration language
    r"\btreated\s+with\b"
    r"|\btreatment\s+with\b"
    r"|\badministration\s+of\b"
    # Reversal by treatment
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

def _induced_by_drug(sent):
    """Return True if the sentence says '<drug> induced' or 'induced by <drug>'."""
    if not re.search(r"\binduced\s+by\b", sent, re.IGNORECASE):
        return False
    return bool(re.search(r"\b(?:" + _DRUG_KW + r")\b", sent, re.IGNORECASE))

def _in_response_to_drug(sent):
    """Return True if 'in response to <drug>'."""
    if not re.search(r"\bin\s+response\s+to\b", sent, re.IGNORECASE):
        return False
    return bool(re.search(r"\b(?:" + _DRUG_KW + r")\b", sent, re.IGNORECASE))

def is_treatment_sentence(sent):
    """Combined check: sentence should be skipped if it describes drug-caused gene changes."""
    if TREATMENT_CAUSAL_RE_NEW.search(sent):
        return True
    if _induced_by_drug(sent):
        return True
    if _in_response_to_drug(sent):
        return True
    return False

GENE_CONTEXT_PATTERNS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bexpress\w*\b",
    r"\bover[\s\-]?express\w*\b",
    r"\bdown[\s\-]?regulat\w*\b",
    r"\bup[\s\-]?regulat\w*\b",
    r"\bdys[\s\-]?regulat\w*\b",
    r"\bregulat\w*\b",
    r"\bgene\s+expression\b",
    r"\bdifferential(ly)?\s+(gene\s+)?expression\b",
    r"\bmrna\s+(expression|level|abundance)\b",
]]

_TIER1_GENES = {"TH","MBP","SST","VIP","RET","TNF","CHAT","CLOCK","SYP","PYY",
                "NPY","CRH","IL6","PER2","SYN1","MAP2","GFAP"}
_BIO_CONTEXT_RE = re.compile(
    r"\b(gene|genes|expression|expressed|mrna|protein|receptor|knockout|"
    r"knockdown|overexpression|transcript|transcription|encoded|encodes|"
    r"regulated|upregulated|downregulated|levels|signaling|signalling|"
    r"pathway|mutation|variant|allele|promoter|locus|immunoreactivity|"
    r"immunostaining|antibody|antibodies|staining|positive|neurons|"
    r"neuronal|intestinal|colonic|hippocampal|cortical|striatal|"
    r"deficiency|deficient|null mutant|heterozygous|homozygous|transgenic|"
    r"knockin|overexpressing|silencing|inhibition|activation)\b", re.IGNORECASE)
_TIER2_ALIASES = {
    "SI":  re.compile(r"\b(sucrase.isomaltase|sucrase|isomaltase|disaccharidase|brush.border enzyme)\b", re.IGNORECASE),
    "GIP": re.compile(r"\b(gastric inhibitory polypeptide|glucose.dependent insulinotropic|incretin|gip receptor|gipr)\b", re.IGNORECASE),
}
BANNED_ACRONYM_CONTEXT = {
    "SI":  [r"suicidal",r"ideation",r"severity",r"international system",r"signal intensity",r"saturation index"],
    "TH":  [r"\bth\s+cell",r"\bth\d+\b",r"\bth\s+hormone",r"thyroid",r"\b\d+th\b"],
}

def _split_sentences(text):
    return re.split(r"(?<=[.!?])\s+|\n", text)

_ABSTRACT_SECTION_RE = re.compile(
    r"(?m)^(?:BACKGROUND|OBJECTIVE[S]?|PURPOSE|INTRODUCTION|AIM[S]?|"
    r"CONTEXT|SUMMARY|METHODS?|RESULTS?|CONCLUSIONS?|ABSTRACT)[:\s]",
    re.IGNORECASE,
)

def _strip_pubmed_header(text):
    """Return only the abstract body, stripping the PubMed citation header."""
    m = _ABSTRACT_SECTION_RE.search(text)
    if m:
        return text[m.start():]
    ai = re.search(r"(?i)\bauthor information\b.*?\n(\n+)", text, re.DOTALL)
    if ai:
        return text[ai.end():]
    return text

def detect_mentioned_genes(text, gene_list):
    tu, tl = text.upper(), text.lower()
    sentences = _split_sentences(text)
    found = set()
    for g in gene_list:
        gu = g.upper()
        pat = re.compile(r"\b" + re.escape(gu) + r"\b")
        if gu in BANNED_ACRONYM_CONTEXT:
            if any(re.search(b, tl) for b in BANNED_ACRONYM_CONTEXT[gu]):
                continue
        if g in _TIER2_ALIASES:
            if pat.search(tu) and _TIER2_ALIASES[g].search(text):
                found.add(g)
        elif g in _TIER1_GENES:
            for sent in sentences:
                if pat.search(sent.upper()) and _BIO_CONTEXT_RE.search(sent):
                    if gu in BANNED_ACRONYM_CONTEXT:
                        if any(re.search(b, sent.lower()) for b in BANNED_ACRONYM_CONTEXT[gu]):
                            continue
                    found.add(g); break
        else:
            if pat.search(tu):
                found.add(g)
    return found

_LINKING_WORDS = {"of","in","was","is","were","are","the","its","their",
                  "a","an","by","for","that","which","to","be","been"}

def detect_gene_context_verbose(text, mentioned_genes, use_new_filter=True):
    """
    Like detect_gene_context but prints a sentence-by-sentence trace.
    Returns True/False (same as production function).
    """
    if _abstract_uses_therapeutic_drug(text):
        print(f"    [ABSTRACT EXCLUDED — therapeutic drug administration detected]")
        return False
    if _abstract_uses_nondrug_intervention(text):
        print(f"    [ABSTRACT EXCLUDED — non-drug intervention detected]")
        return False
    if _abstract_attributes_expression_to_treatment(text):
        print(f"    [ABSTRACT EXCLUDED — expression attributed to treatment (generic)]")
        return False
    if _abstract_uses_epigenetics(text):
        print(f"    [ABSTRACT EXCLUDED — epigenetics/miRNA study]")
        return False
    if _abstract_is_behavioral_primary(text):
        print(f"    [ABSTRACT EXCLUDED — behavioral primary, gene expression secondary]")
        return False
    if _is_genetic_association_study(text):
        print(f"    [ABSTRACT EXCLUDED — genetic association/GWAS study, no expression measurement]")
        return False
    if not _has_disease_model(text):
        print(f"    [ABSTRACT EXCLUDED — no disease model detected]")
        return False
    if not _EXPRESSION_MEASUREMENT_RE.search(_get_results_section(text)):
        print(f"    [ABSTRACT EXCLUDED — no gene expression measurement signal in results section]")
        return False

    gene_terms = [r"genes?"] + [re.escape(g) for g in mentioned_genes]
    gene_any_pat = re.compile(r"\b(?:" + "|".join(gene_terms) + r")\b", re.IGNORECASE)

    result = False
    for sent in _split_sentences(text):
        sent = sent.strip()
        if not sent:
            continue
        has_ctx = any(p.search(sent) for p in GENE_CONTEXT_PATTERNS)
        if not has_ctx:
            continue  # no regulation word → skip silently

        old_excluded = bool(TREATMENT_CAUSAL_RE_OLD.search(sent))
        new_excluded = is_treatment_sentence(sent)

        # Check proximity
        words = sent.split()
        n = len(words)
        word_starts = []
        cursor = 0
        for w in words:
            idx = sent.find(w, cursor)
            word_starts.append(idx)
            cursor = idx + len(w)

        passes_proximity = False
        for pat in GENE_CONTEXT_PATTERNS:
            for m in pat.finditer(sent):
                reg_words = [wi for wi, ws in enumerate(word_starts)
                             if ws <= m.start() < ws + len(words[wi])
                             or m.start() <= ws < m.end()]
                if not reg_words:
                    continue
                reg_lo, reg_hi = min(reg_words), max(reg_words)
                for j in (reg_lo - 1, reg_hi + 1):
                    if 0 <= j < n and gene_any_pat.search(words[j]):
                        passes_proximity = True
                for gene_j, mid_j in ((reg_lo-2, reg_lo-1),(reg_hi+2, reg_hi+1)):
                    if 0 <= gene_j < n and 0 <= mid_j < n:
                        mid = words[mid_j].lower().rstrip(".,;:")
                        if mid in _LINKING_WORDS and gene_any_pat.search(words[gene_j]):
                            passes_proximity = True
                for gene_j, m1j, m2j in ((reg_lo-3,reg_lo-2,reg_lo-1),(reg_hi+3,reg_hi+1,reg_hi+2)):
                    if all(0 <= j < n for j in (gene_j, m1j, m2j)):
                        m1 = words[m1j].lower().rstrip(".,;:")
                        m2 = words[m2j].lower().rstrip(".,;:")
                        if (m1 in _LINKING_WORDS or m2 in _LINKING_WORDS) and gene_any_pat.search(words[gene_j]):
                            passes_proximity = True

        excluded = new_excluded if use_new_filter else old_excluded
        diverges = (old_excluded != new_excluded)

        # Determine the sentence outcome under each filter
        old_result = passes_proximity and not old_excluded
        new_result = passes_proximity and not new_excluded

        # Only print sentences where something interesting happens
        interesting = passes_proximity or diverges
        if interesting:
            tag_old = "EXCLUDED(old)" if old_excluded else "PASS(old)"
            tag_new = "EXCLUDED(new)" if new_excluded else "PASS(new)"
            flag = " ◄ FILTER DIVERGES" if diverges else ""
            print(f"    [{tag_old} | {tag_new}]{flag}")
            print(f"    Proximity: {'YES' if passes_proximity else 'NO'}")
            # Truncate long sentences
            preview = sent[:200] + ("…" if len(sent) > 200 else "")
            print(f"    Sentence: {preview}")
            print()

        if new_result:
            result = True

    return result


def fetch_small_batch(query, year, n):
    fq = (f"{query} AND animals[MeSH Terms] AND {year}[dp] "
          f"NOT (clinical trial[pt] OR randomized controlled trial[pt] "
          f"OR meta-analysis[pt] OR review[pt])")
    handle = Entrez.esearch(db="pubmed", term=fq, retmax=n)
    record = Entrez.read(handle); handle.close()
    ids = record["IdList"][:n]
    if not ids:
        print("No results."); return []
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
    data = handle.read(); handle.close()
    abstracts = re.split(r"\n\n(?=\d+\.\s)", data)
    return [a for a in abstracts if len(a.strip()) > 50]


if __name__ == "__main__":
    print(f"\n{'='*70}")
    print(f"  DEBUG RUN — {DEBUG_DISORDER} — {DEBUG_YEAR} — {DEBUG_N} abstracts")
    print(f"{'='*70}\n")

    abstracts = fetch_small_batch(DEBUG_QUERY, DEBUG_YEAR, DEBUG_N)
    print(f"Fetched {len(abstracts)} abstracts.\n")

    # ── Optional LLM-gate smoke test ─────────────────────────────────────────
    # USE_LLM_GATE=1 bypasses the verbose sentence trace and instead batch-qualifies
    # the abstracts through llm_qualify (the same gate the production scanner uses),
    # so you can eyeball the LLM's qualifies/reason/tissue_hint on a tiny fetch.
    if os.environ.get("USE_LLM_GATE", "0") == "1":
        try:
            from llm_qualify import qualify_batch_sync
        except Exception as e:
            print(f"llm_qualify unavailable ({e}); set USE_LLM_GATE=0 to use the trace.")
            raise SystemExit
        items = []
        for i, abstract in enumerate(abstracts):
            body = _strip_pubmed_header(abstract)
            mentioned = detect_mentioned_genes(body, ALL_CANDIDATE_GENES)
            pmid_m = re.search(r'\bPMID:\s*(\d+)', abstract)
            iid = (pmid_m.group(1) if pmid_m
                   else f"dbg{i}:" + hashlib.sha256(body.encode()).hexdigest()[:16])
            items.append((iid, body, sorted(mentioned)))
        judgments = qualify_batch_sync(items)
        print(f"{'─'*70}\nLLM gate — {len(judgments)} judgments:\n{'─'*70}")
        for i, (iid, _body, _mg) in enumerate(items):
            title = (abstracts[i].strip().split("\n")[1]
                     if len(abstracts[i].strip().split("\n")) > 1 else "")
            j = judgments.get(iid)
            if j:
                print(f"#{i+1} qualifies={j.qualifies}  tissue_hint={j.tissue_hint}\n"
                      f"   reason: {j.reason}\n   title: {title.strip()}\n")
            else:
                print(f"#{i+1} (no judgment)\n")
        raise SystemExit

    for i, abstract in enumerate(abstracts):
        title_line = abstract.strip().split("\n")[1] if len(abstract.strip().split("\n")) > 1 else ""
        print(f"{'─'*70}")
        print(f"Abstract #{i+1}: {title_line.strip()}")

        body = _strip_pubmed_header(abstract)
        mentioned = detect_mentioned_genes(body, ALL_CANDIDATE_GENES)
        print(f"  Genes detected: {sorted(mentioned) or '(none)'}")

        # Check for drug-related words in the abstract body (red flag)
        drug_in_abstract = bool(re.search(
            r"\b(?:drug[s]?|medication[s]?|antidepressant[s]?|antipsychotic[s]?|"
            r"haloperidol|clozapine|risperidone|fluoxetine|phencyclidine|pcp|"
            r"amphetamine|methylphenidate|ketamine|lithium|valproate|dexamethasone|"
            r"treated\s+with|treatment\s+with|administered)\b",
            body, re.IGNORECASE))
        if drug_in_abstract:
            print(f"  ⚠  Drug/treatment keywords found in abstract body")

        if not mentioned:
            print(f"  → Skipped (no candidate genes found)\n")
            continue

        # Print the abstract body so we can spot parsing issues
        print(f"\n  [Abstract body preview]")
        for line in body.strip().split("\n")[:12]:
            print(f"    {line}")
        print()

        print(f"\n  Sentence analysis (sentences with gene-context words only):")
        result = detect_gene_context_verbose(body, mentioned, use_new_filter=True)
        print(f"  ► gene_context (NEW filter) = {result}")

        # Compute old-filter result with full proximity logic (1/2/3 words)
        old_result = False
        gene_terms = [r"genes?"] + [re.escape(g) for g in mentioned]
        gene_any_pat = re.compile(r"\b(?:" + "|".join(gene_terms) + r")\b", re.IGNORECASE)
        for sent in _split_sentences(body):
            if not any(p.search(sent) for p in GENE_CONTEXT_PATTERNS):
                continue
            if TREATMENT_CAUSAL_RE_OLD.search(sent):
                continue
            words = sent.split(); n2 = len(words)
            word_starts = []
            cursor = 0
            for w in words:
                idx = sent.find(w, cursor); word_starts.append(idx); cursor = idx + len(w)
            for pat in GENE_CONTEXT_PATTERNS:
                for m in pat.finditer(sent):
                    reg_words = [wi for wi, ws in enumerate(word_starts)
                                 if ws <= m.start() < ws + len(words[wi]) or m.start() <= ws < m.end()]
                    if not reg_words: continue
                    reg_lo, reg_hi = min(reg_words), max(reg_words)
                    for j in (reg_lo-1, reg_hi+1):
                        if 0 <= j < n2 and gene_any_pat.search(words[j]):
                            old_result = True
                    for gene_j, mid_j in ((reg_lo-2,reg_lo-1),(reg_hi+2,reg_hi+1)):
                        if 0 <= gene_j < n2 and 0 <= mid_j < n2:
                            mid = words[mid_j].lower().rstrip(".,;:")
                            if mid in _LINKING_WORDS and gene_any_pat.search(words[gene_j]):
                                old_result = True
                    for gene_j, m1j, m2j in ((reg_lo-3,reg_lo-2,reg_lo-1),(reg_hi+3,reg_hi+1,reg_hi+2)):
                        if all(0 <= j < n2 for j in (gene_j, m1j, m2j)):
                            m1 = words[m1j].lower().rstrip(".,;:")
                            m2 = words[m2j].lower().rstrip(".,;:")
                            if (m1 in _LINKING_WORDS or m2 in _LINKING_WORDS) and gene_any_pat.search(words[gene_j]):
                                old_result = True
        if result != old_result:
            print(f"  ◄ FILTER CHANGED THE OUTCOME: old={old_result} → new={result}")
        print()
