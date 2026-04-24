# LCP1-NET: Network-Level Characterisation of L-Plastin Co-expression Dysregulation in Colorectal Cancer

<div align="center">

![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.3%2B-276DC3?style=flat-square&logo=r&logoColor=white)
![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-EE4C2C?style=flat-square&logo=pytorch&logoColor=white)
![Neo4j](https://img.shields.io/badge/Neo4j-5.x-008CC1?style=flat-square&logo=neo4j&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-0F6E56?style=flat-square)
![Status](https://img.shields.io/badge/Status-Active%20Research-orange?style=flat-square)

**Masters Thesis Project · University of Padova × Dublin City University**

*An unsupervised machine learning pipeline that detects gene co-expression network breakdown in colorectal cancer — with a specific focus on the oncogenic role of L-Plastin (LCP1)*

[Overview](#overview) · [Architecture](#architecture) · [Installation](#installation) · [Data](#data) · [Usage](#usage) · [Results](#results) · [Citation](#citation)

---

</div>

## Overview

Traditional RNA-seq analysis asks: *"Is this gene up or down?"*

This project asks something fundamentally different: *"Has this entire gene network lost coordination — and is LCP1 driving that breakdown?"*

**LCP1-NET** is a computational pipeline that:

1. Constructs a normal colorectal tissue gene co-expression network from GTEx and TCGA healthy samples
2. Identifies which module L-Plastin (LCP1) belongs to and characterises its co-expression neighbourhood
3. Trains an autoencoder on normal network patterns to learn what healthy gene coordination looks like
4. Scores individual patient tumours (TCGA-COAD and TCGA-READ) for network-level anomalies — per patient, per module
5. Connects anomaly scores to biological pathways through a Neo4j graph database
6. Tests whether LCP1 network dysregulation is correlated with expression level, tumour stage, and patient survival

> **Why LCP1?** L-Plastin is an actin-bundling protein normally restricted to haematopoietic cells. Its aberrant expression in colorectal cancer drives invasion, E-cadherin loss, and poor prognosis — but its network-level behaviour across a large patient cohort has never been characterised computationally.

---

## Key Findings

| Finding | Result |
|---|---|
| LCP1 module | `ME[INSERT]` · `[N]` member genes |
| LCP1 hub gene status | kME = `[INSERT]` (hub: `[yes/no]`) |
| COAD anomalous patients | `[N]`/`[N]` (`[N]`%) at 2σ threshold |
| READ anomalous patients | `[N]`/`[N]` (`[N]`%) at 2σ threshold |
| LCP1 expression vs anomaly | r = `[INSERT]`, p = `[INSERT]` |
| Stage dependence | Spearman ρ = `[INSERT]`, p = `[INSERT]` |
| Overall survival (combined) | log-rank p = `[INSERT]` |
| Top enriched pathway | `[INSERT]` |

> Replace bracketed values with your computed results after running the full pipeline.

---

## Architecture

The pipeline runs in three sequential stages. Each stage feeds into the next.

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1 — R Pipeline (WGCNA)                                   │
│                                                                  │
│  GTEx Sigmoid Colon  ──┐                                        │
│  GTEx Transverse Colon ┼──► ComBat-seq ──► WGCNA ──► Modules   │
│  TCGA-COAD normals   ──┤                       │                │
│  TCGA-READ normals   ──┘                       ▼                │
│                                         Gene-module map         │
│                                         Normal eigengenes       │
│                                         Hub gene rankings       │
└──────────────────────────────────┬──────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2 — Python ML Pipeline                                   │
│                                                                  │
│  Normal eigengenes ──► Autoencoder training ──► Thresholds      │
│                              (normal only)                       │
│  TCGA-COAD tumours ──► Project eigengenes ──► Anomaly scores    │
│  TCGA-READ tumours ──► Project eigengenes ──► Anomaly scores    │
│                                                                  │
│  LCP1-specific: expression correlation · stage dependence       │
│                 EMT co-dysregulation · survival analysis        │
└──────────────────────────────────┬──────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3 — Neo4j Graph Database                                 │
│                                                                  │
│  Gene nodes ──► Module nodes ──► Pathway nodes ──► GO terms     │
│                                                                  │
│  Patient nodes ──► SCORED_MODULE ──► Module nodes               │
│                                                                  │
│  Cypher queries: LCP1 neighbourhood · pathway enrichment        │
│                  CRC driver validation · patient reports        │
└─────────────────────────────────────────────────────────────────┘
```

### Why this combination?

| Tool | Role | Why not an alternative? |
|---|---|---|
| **WGCNA** | Co-expression network construction | Data-driven modules capture actual tissue biology; fixed gene sets (GSVA) miss novel relationships |
| **Autoencoder** | Anomaly detection | Unsupervised — no disease labels needed; learns normal coordination patterns not individual expression levels |
| **Neo4j** | Biological knowledge integration | Native graph queries over gene-pathway-patient relationships; SQL would require complex joins for every biological question |
| **recount3** | Uniform data acquisition | All sources (GTEx + TCGA) processed through one pipeline — minimises technical batch effects |

---

## Repository Structure

```
LCP1-NET/
│
├── r_pipeline/
│   └── complete_pipeline.R          # Full WGCNA pipeline — data to modules
│
├── ml_pipeline/
│   ├── 00_config.py                 # All hyperparameters and paths
│   ├── 01_preprocess.py             # Load WGCNA outputs, project tumour eigengenes
│   ├── 02_train_autoencoder.py      # Train and validate autoencoder
│   ├── 03_score_anomalies.py        # Score COAD and READ tumour samples
│   ├── 04_lcp1_analysis.py          # LCP1-specific biological analyses
│   ├── 05_clinical_validation.py    # Kaplan-Meier survival analysis
│   ├── run_all.py                   # Master script — runs all steps
│   └── models/
│       └── autoencoder.py           # CRCAutoencoder class definition
│
├── neo4j_pipeline/
│   ├── 00_neo4j_config.py           # Connection settings and gene identifiers
│   ├── 01_build_graph.py            # Load gene and module nodes from WGCNA
│   ├── 02_fetch_biological_data.py  # Download KEGG and GO annotations
│   ├── 03_load_anomaly_scores.py    # Load ML scores as Patient nodes
│   ├── 04_lcp1_queries.py           # 12 Cypher queries for biological results
│   └── run_neo4j_pipeline.py        # Master script
│
├── data/
│   ├── raw/                         # Downloaded counts (gitignored — large)
│   ├── wgcna/                       # WGCNA outputs from R pipeline
│   └── processed/                   # Tumour expression matrices
│
├── ml_outputs/
│   ├── models/                      # Trained autoencoder weights + scaler
│   ├── scores/                      # Per-patient anomaly score CSVs
│   └── figures/                     # All publication-quality plots
│
├── neo4j_results/                   # CSV outputs from Cypher queries
│
├── requirements.txt
├── environment.R
└── README.md
```

---

## Installation

### Prerequisites

- R ≥ 4.3.0
- Python ≥ 3.9
- Neo4j Desktop ≥ 5.x (free for research)
- 16 GB RAM recommended
- ~20 GB free disk space for data

### 1. Clone the repository

```bash
git clone https://github.com/[YOUR_USERNAME]/LCP1-NET.git
cd LCP1-NET
```

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

**requirements.txt:**
```
torch>=2.0.0
numpy>=1.24.0
pandas>=2.0.0
scikit-learn>=1.3.0
scipy>=1.11.0
matplotlib>=3.7.0
seaborn>=0.12.0
lifelines>=0.27.0
neo4j>=5.0.0
requests>=2.31.0
tqdm>=4.65.0
```

### 3. Install R dependencies

```r
source("environment.R")
```

**environment.R:**
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "recount3", "DESeq2", "sva",
    "WGCNA", "SummarizedExperiment", "biomaRt"
))

install.packages(c("flashClust", "ggplot2", "gridExtra"))
```

### 4. Set up Neo4j

1. Download [Neo4j Desktop](https://neo4j.com/download/) — free for research
2. Create a new local database
3. Set password to `crc_pipeline` (or update `neo4j_pipeline/00_neo4j_config.py`)
4. Start the database before running Stage 3

---

## Data

All data is publicly available. No special access is required beyond standard registration.

### Normal reference (Stage 1)

| Source | Tissue | Access | Approx. samples |
|---|---|---|---|
| GTEx V10 | Colon — Sigmoid | [gtexportal.org](https://gtexportal.org) — free registration | ~320 |
| GTEx V10 | Colon — Transverse | [gtexportal.org](https://gtexportal.org) — free registration | ~280 |
| TCGA-COAD | Matched normal tissue | Via recount3 — no application needed | ~40 |
| TCGA-READ | Matched normal tissue | Via recount3 — no application needed | ~10 |

> GTEx requires a Data Use Agreement. Apply at gtexportal.org — approval typically takes 24–48 hours.

### Tumour validation (Stage 2)

| Cohort | Description | Access |
|---|---|---|
| TCGA-COAD | Colon adenocarcinoma primary tumours | Via recount3 |
| TCGA-READ | Rectal adenocarcinoma primary tumours | Via recount3 |

All data is downloaded programmatically by the R pipeline. No manual file management required.

---

## Usage

### Stage 1 — WGCNA Network Construction (R)

```r
source("r_pipeline/complete_pipeline.R")
# Runtime: 2-4 hours
# Outputs: data/wgcna/*.csv
```

Check LCP1 immediately after Stage 1:
```r
lcp1_module <- mergedColors[colnames(expression.data) == "ENSG00000159111"]
cat("LCP1 module:", lcp1_module, "\n")
```

### Stage 2 — Machine Learning Pipeline (Python)

```bash
cd ml_pipeline
python run_all.py
# Runtime: 20-60 minutes total
# Outputs: ml_outputs/
```

Or run individual steps:
```bash
python 01_preprocess.py          # ~5 min — project tumour eigengenes
python 02_train_autoencoder.py   # ~10-30 min — train on normal data only
python 03_score_anomalies.py     # ~5 min — score COAD and READ patients
python 04_lcp1_analysis.py       # ~5 min — LCP1-specific analyses
python 05_clinical_validation.py # ~5 min — Kaplan-Meier survival
```

### Stage 3 — Neo4j Graph Database (Python)

```bash
# Ensure Neo4j Desktop is running first
cd neo4j_pipeline
python run_neo4j_pipeline.py
# Runtime: 30-90 minutes (API calls are the bottleneck)
# Outputs: neo4j_results/Q1-Q12 CSV files
```

View the LCP1 network neighbourhood visually at `http://localhost:7474`:
```cypher
MATCH (lcp1:LCP1Gene)-[:BELONGS_TO]->(m:Module)
MATCH (partner:Gene)-[:BELONGS_TO]->(m)
WHERE partner.is_hub = true
OPTIONAL MATCH (partner)-[:IN_PATHWAY]->(p:Pathway)
WHERE p.is_lcp1_associated = true
RETURN lcp1, m, partner, p
LIMIT 50
```

---

## Autoencoder Architecture

```
Input        n_modules (e.g. 18 eigengenes)
    │
Encoder      128 → 64 → 8 neurons   BatchNorm · ReLU · Dropout(0.2)
    │
Latent        8 neurons · Tanh       compressed "normal" representation
    │
Decoder       64 → 128 → n_modules  BatchNorm · ReLU · Dropout(0.2)

Training loss:  MSE(input, reconstruction)
Anomaly score:  mean((input − reconstruction)²) per sample
Threshold:      μ + 2σ of normal test set reconstruction errors
```

The model is trained **exclusively on normal tissue**. It never sees tumour data during training. High reconstruction error on a new sample means that sample's gene network pattern was not learned as normal — that is the anomaly signal.

---

## Graph Database Schema

```cypher
(:Gene {gene_id, gene_symbol, kME, is_hub, is_lcp1,
        is_emt_marker, is_crc_driver})
    -[:BELONGS_TO]->
(:Module {name, colour, n_genes, is_lcp1_module})

(:Gene)-[:IN_PATHWAY]->(:Pathway {pathway_id, name, is_lcp1_associated})

(:Gene)-[:HAS_GO_TERM {evidence}]->(:GOTerm {go_id, name, aspect})

(:Patient {patient_id, cohort, total_anomaly_score,
           lcp1_module_score, tumor_stage, vital_status})
    -[:SCORED_MODULE {anomaly_score, is_anomalous}]->
(:Module)
```

---

## Comparison to Existing Tools

| Tool | What it detects | Gap this project fills |
|---|---|---|
| **OUTRIDER** (Brechtmann et al., 2018) | Single-gene expression outliers; treats co-variation as a confound | We treat co-variation as the signal — detecting coordination breakdown not expression deviation |
| **GSVA / ssGSEA** | Gene set activity using fixed curated sets | We use data-driven co-expression modules from the specific tissue studied |
| **multiWGCNA** | Differential co-expression between predefined groups | We score individual patient samples against a pre-trained normal reference — enabling prospective analysis |

---

## Troubleshooting

<details>
<summary>LCP1 not found in gene-module assignments</summary>

LCP1 may have been filtered out during low-expression gene removal. Lower the threshold:
```r
min_samples <- ceiling(ncol(normal_combined) * 0.10)  # was 0.20
```

</details>

<details>
<summary>Batch correction does not converge</summary>

Try preserving tissue biology when correcting:
```r
corrected_counts <- ComBat_seq(
    counts = normal_filtered,
    batch  = batch_labels,
    group  = tissue_labels   # Preserve colon vs rectum biology
)
```

</details>

<details>
<summary>Neo4j connection refused</summary>

Ensure the database is started in Neo4j Desktop (green status indicator). Default bolt port is `7687`. Test with:
```python
from neo4j import GraphDatabase
driver = GraphDatabase.driver("bolt://localhost:7687",
                               auth=("neo4j", "crc_pipeline"))
with driver.session() as s:
    print(s.run("RETURN 1").single()[0])
```

</details>

<details>
<summary>Autoencoder overfitting (validation loss diverges)</summary>

Reduce model complexity in `ml_pipeline/00_config.py`:
```python
HIDDEN_1     = 64     # was 128
HIDDEN_2     = 32     # was 64
LATENT_DIM   = 4      # was 8
DROPOUT_RATE = 0.3    # was 0.2
```

</details>

---

## Extending the Pipeline

**Different gene of interest:**
```python
# ml_pipeline/00_config.py
LCP1_ENSEMBL = "ENSG00000XXXXXX"
LCP1_SYMBOL  = "YOURGENE"
```

**Different cancer type:**
```r
# r_pipeline/complete_pipeline.R
coad_info <- subset(tcga_projects, project == "BRCA")  # Breast cancer
read_info  <- subset(tcga_projects, project == "LUAD") # Lung adenocarcinoma
```

---

## Citation

```bibtex
@mastersthesis{lcp1net2025,
  title  = {Network-Level Characterisation of L-Plastin (LCP1)
            Co-expression Dysregulation in Colorectal Cancer},
  author = {[YOUR NAME]},
  school = {University of Padova / Dublin City University},
  year   = {2025},
  note   = {Quantitative and Computational Biosciences}
}
```

### Dependencies to cite

| Tool | Reference |
|---|---|
| WGCNA | Langfelder & Horvath (2008). *BMC Bioinformatics* |
| recount3 | Wilks et al. (2021). *Genome Biology* |
| DESeq2 | Love et al. (2014). *Genome Biology* |
| ComBat-seq | Zhang et al. (2020). *NAR Genomics and Bioinformatics* |
| GTEx | GTEx Consortium (2020). *Science* |
| TCGA | Cancer Genome Atlas Research Network (2012). *Nature* |
| PyTorch | Paszke et al. (2019). *NeurIPS* |

---

## Acknowledgements

- **Dublin City University** — research facilities and supervision
- **University of Padova** — academic programme in Quantitative and Computational Biosciences
- **GTEx Consortium** — for making normal tissue RNA-seq data publicly available
- **TCGA Research Network** — for making cancer genomics data publicly available
- **recount3 team** — for uniformly processing hundreds of thousands of RNA-seq samples

---

## License

MIT License — see [LICENSE](LICENSE) for details.

Data from GTEx and TCGA is used in accordance with their respective data use agreements. Biological annotations from KEGG and Gene Ontology are used for non-commercial academic research.

---

<div align="center">

**Masters Thesis · Quantitative and Computational Biosciences**
**University of Padova × Dublin City University · 2025**

*Built to advance our understanding of L-Plastin's network-level role in colorectal cancer*

</div>
