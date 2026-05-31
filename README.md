# 🔬 Biomarker Quantification in Breast Cancer — Xenium In Situ

**Reproducing Figures 1A, 1B, 1C, 1E, 1F, 2C, 2F, 2G from:**  
> Janesick A, Kravitz SN, Stauffer W, Valencia M, Taylor SEB, et al.  
> *Biomarker Quantification in Breast Cancer using Xenium In Situ.*  
> bioRxiv (2025). https://doi.org/10.64898/2025.12.08.692193

**Dataset:** [10x Genomics — Xenium FFPE Human Breast Biomarkers](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers)  
**Platform:** Google Colab (Python 3, GPU not required)

---

## 📋 Table of Contents

1. [Understanding the Paper](#1-understanding-the-paper)  
   - [Purpose](#11-purpose)  
   - [Background & Existing Work](#12-background--existing-work)  
   - [Methodology](#13-detailed-description-of-methodology)  
   - [Results](#14-results)  
2. [Project Structure](#2-project-structure)  
3. [Input Data Files](#3-input-data-files)  
4. [Setup & Installation](#4-setup--installation)  
5. [Step-by-Step Code Walkthrough](#5-step-by-step-code-walkthrough)  
6. [Output Files & Figures](#6-output-files--figures)  
7. [Key Tables & Data Summaries](#7-key-tables--data-summaries)  
8. [Troubleshooting](#8-troubleshooting)  
9. [Dependencies](#9-dependencies)

---

## 1. Understanding the Paper

### 1.1 Purpose

Cancer diagnosis has a big problem: you can detect gene expression, but how do you know whether a gene is *actually* behaving abnormally, or whether the measurement itself is just noisy? This is the core question Janesick et al. (2025) set out to answer. The paper has two central goals:

1. **Find reliable "housekeeping" (HK) genes** — genes whose expression stays stable regardless of how aggressive the cancer is, so they can be used as internal rulers to normalize other measurements.
2. **Identify biomarkers that actually track tumor aggressiveness** — genes whose expression goes up or down predictably as Ductal Carcinoma In Situ (DCIS) progresses from low-grade to high-grade, and potentially toward invasion.

The technology they use — **Xenium In Situ** — makes this possible at single-cell resolution, inside intact FFPE (formalin-fixed, paraffin-embedded) breast tissue. Rather than homogenizing tissue and measuring average gene expression from millions of cells mixed together, Xenium reads individual transcripts (RNA molecules) at their precise spatial location in the tissue. This means you can look at tumor cells, myoepithelial cells, and stromal fibroblasts separately — and even map where genes are expressed relative to duct boundaries, down to a 30 µm scale.

---

### 1.2 Background & Existing Work

#### Why DCIS matters
DCIS is a pre-invasive breast cancer — malignant cells are confined within the milk ducts and have not yet broken through into surrounding tissue. Roughly 20–50% of untreated DCIS cases eventually become invasive, but we currently cannot reliably predict which ones will. If we could, many patients could be spared from aggressive treatment. This makes accurate molecular biomarkers critically important.

#### The normalization problem in spatial transcriptomics
Any RNA measurement platform — whether it's bulk RNA-seq, RT-qPCR, or in situ — produces raw counts that are affected by technical noise: tissue preparation quality, FFPE degradation, RNA capture efficiency, and sample-to-sample variation. To compare gene expression meaningfully across different tissue sections or patients, you need a normalization strategy.

The most common approach, borrowed from clinical diagnostic assays like **Oncotype Dx** (a breast cancer recurrence test), is to divide target gene counts by the counts of a set of "housekeeping" genes that are assumed to be stably expressed across all samples. Oncotype Dx uses five HK genes: *ACTB*, *GAPDH*, *GUSB*, *RPLP0*, and *TFRC*.

#### The problem with classic HK genes at single-cell resolution
Traditional HK gene selection was done using bulk RNA measurements from mixed cell populations. But at single-cell resolution — which Xenium provides — you can now see something that was previously invisible: **some of those "stable" HK genes are not actually stable at all when you look within specific cell types**.

For example, *GAPDH* is a glycolytic enzyme. As tumors become more aggressive, they ramp up glycolysis (the Warburg effect), which means *GAPDH* expression *itself* increases in more aggressive tumor cells. Using *GAPDH* as a normalization reference would therefore mask or even reverse the apparent expression of tumor biomarkers. Similarly, *TFRC* and *GUSB* showed significant variability across different tumor subtypes in this study.

#### What spatial in situ adds
Earlier work by the same group (Janesick et al., Nature Communications 2023) demonstrated the power of combining Xenium with Visium spatial transcriptomics and single-cell RNA-seq on the same tissue block. That paper characterized cell type composition of the breast tumor microenvironment at high resolution. The 2025 paper extends this to focus specifically on biomarker quantification and normalization — a more practically translatable problem for clinical diagnostics.

---

### 1.3 Detailed Description of Methodology

The study analyzed **12 FFPE human breast tissue samples** (8 DCIS, 3 invasive ductal carcinoma, 1 normal) using a **custom 280-gene Xenium panel** that included:
- Breast tissue cell type markers (epithelial, myoepithelial, stromal, immune)
- 45 candidate housekeeping genes from literature and internal databases
- Known clinical biomarkers (ER/ESR1, PR/PGR, HER2/ERBB2, MKI67)
- Basement membrane (BM) genes relevant to invasion potential

#### Step 1: Cell segmentation and cell typing
Each tissue section was processed through Xenium's on-instrument pipeline (v4.0.0), which produces cell boundaries, nucleus centroids, and per-cell transcript counts. Cells were filtered to retain those with >10 total transcripts and >5 unique genes detected. Normalization used SCTransform, followed by PCA, neighbor graph construction, UMAP, and Leiden clustering (resolution 0.3). Cell types were annotated using known marker genes.

#### Step 2: Assessing HK gene stability using Coefficient of Variation (CV)
The CV (standard deviation / mean) is a measure of relative variability. A low CV means a gene's expression is consistent from cell to cell within a given cell type across all 12 sections. A high CV means it varies a lot.

The authors computed **intra-section CV** (within each tissue section, across cells of the same type) and **inter-section CV** (across all 12 sections, comparing cluster-level means). This two-level analysis was critical: a gene might be stable within one section but vary dramatically between sections (making it useless for cross-sample normalization).

Through this analysis, four genes emerged as the best HK genes for tumor cells:
- **EEF1G** — translation elongation factor, low CV across all cell types
- **EEF2** — translation elongation factor, similarly stable
- **MALAT1** — a highly abundant long non-coding RNA, stably expressed
- **RPLP0** — ribosomal protein, stable and well-supported in literature

#### Step 3: Choosing which biomarkers to normalize
The study then identified genes with **high inter-section CV** — meaning they vary significantly across different sections — as candidate biomarkers. These are the genes that *should* vary because they track tumor biology. Crucially, the normalization strategy must preserve this variance; a bad HK gene set would flatten the differences between tumor grades.

Key biomarkers identified:
- **LDHA** and **SDC1**: increase with tumor grade (linked to the Warburg effect and aggressive behavior)
- **SFRP1** and **PIGR**: decrease with tumor grade (tumor suppressors)

#### Step 4: Validating normalization — HK gene comparison
The authors directly tested whether normalizing LDHA by the Oncotype Dx HK panel (including GAPDH and GUSB) vs. their refined HK panel (EEF1G, EEF2, MALAT1, RPLP0) made a difference. The answer was stark: GAPDH normalization actually reversed the apparent LDHA trend, making high-grade DCIS look similar to benign tissue. This is because GAPDH increases with tumor aggressiveness (Warburg effect), partially canceling out the LDHA signal.

#### Step 5: Spatial biomarker analysis — myoepithelial cells
The study then zoomed into myoepithelial cells — the ring of contractile cells surrounding the milk duct that forms a physical barrier against invasion. They compared myoepithelial cells adjacent to **normal ducts** versus those adjacent to **DCIS ducts**, looking at expression of basement membrane (BM) genes.

Key finding: **LAMC2** (Laminin subunit gamma 2) is significantly upregulated in tumor-associated myoepithelial cells compared to normal ones, while other myoepithelial markers (KRT14, COL17A1) remain relatively stable. LAMC2 protein is cleaved proteolytically in a tumorigenic environment, weakening cell-matrix adhesion and potentially enabling invasion. This makes LAMC2 an early molecular indicator of invasive potential — potentially detectable before histological changes are visible.

#### Step 6: Cell-agnostic spatial analysis of peritumoral MMP11
The final methodological innovation addresses a fundamental limitation of cell-based spatial analysis: **not all transcripts get assigned to cells**. In the dense extracellular matrix (ECM) surrounding ducts, stromal cells (fibroblasts, endothelial cells) have irregular, thin shapes that are hard to segment. Transcripts in the ECM space between cells are simply "unassigned" — and this creates an underestimation of ECM-related genes.

To overcome this, the authors used a **cell-agnostic approach** for *MMP11*, a matrix metalloproteinase expressed in the peritumoral stroma:
1. Manually drew 16 ROI (Region of Interest) polygons around ducts in Xenium Explorer (8 around high-MMP11 ducts, 8 around low-MMP11 ducts) across two tissue sections
2. Built a **30 µm periductal ring** outside each duct boundary
3. Counted all MMP11 transcripts within the ring — regardless of cell assignment
4. Expressed counts as **transcripts per µm²** to normalize for ring area

This approach captures the full MMP11 signal including unassigned transcripts, making it more sensitive than per-cell quantification alone.

---

### 1.4 Results

| Finding | Key Gene(s) | Significance |
|---|---|---|
| Best HK genes for tumor cells | EEF1G, EEF2, MALAT1, RPLP0 | Low CV intra- and inter-section; suitable for normalization |
| Oncotype Dx HK genes problematic | TFRC, GUSB, GAPDH | High CV; GAPDH is Warburg-affected and can reverse biomarker trends |
| Grade-increasing biomarker | LDHA, SDC1 | Increase with DCIS progression; LDHA is Warburg-associated |
| Grade-decreasing biomarker | SFRP1, PIGR | Decrease with DCIS grade; both are tumor suppressors |
| BM invasion marker | LAMC2 | Upregulated in myoepithelial cells around DCIS vs. normal ducts; early invasion signal |
| Peritumoral invasion signal | MMP11 | High periductal density associated with proliferating tumor cells (MKI67+, CCNB1+) and aggressive gene signatures (CEACAM6, S100P, NAT1) |

The most clinically important takeaway is that **selecting appropriate HK genes at single-cell resolution matters enormously**. Using the wrong ones can actively mislead the analysis — suppressing real biological signals or reversing apparent trends. The paper provides a rigorous framework for doing this systematically on any Xenium dataset.

---

## 2. Project Structure

```
xenium_breast_cancer/
│
├── Xenium_Breast_Cancer_Biomarker_Analysis.ipynb   # Main analysis notebook
│
├── xenium_data/                    # Input data (Google Drive)
│   ├── cell_feature_matrix.zarr.zip   (25.1 MB)
│   ├── cells.zarr.zip                 (206.2 MB)
│   ├── transcripts.zarr.zip           (1.18 GB)
│   └── gene_panel.json                (153 KB)
│
└── xenium_outputs/                 # Generated by notebook
    ├── qc_distributions.png
    ├── fig1A_hk_cv_violin.png
    ├── fig1B_deg_cv_violin.png
    ├── fig1C_cv_scatter.png
    ├── fig1E_LDHA_raw.png
    ├── fig1F_LDHA_normalized.png
    ├── fig2C_myoepithelial_markers.png
    ├── fig2F_MMP11_raw_density.png
    ├── fig2G_MMP11_normalized_density.png
    ├── xenium_annotated.h5ad
    ├── hk_gene_cv_table.csv
    └── tableS3_stable_genes.csv
```

---

## 3. Input Data Files

The dataset is downloaded from the [10x Genomics dataset page](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers), specifically from the **S1-Top** tissue section:

```
https://cf.10xgenomics.com/samples/xenium/4.0.0/Human_Breast_Biomarkers_S1_Top/
Human_Breast_Biomarkers_S1_Top_xe_outs.zip
```

After extracting, the following four files are used directly:

| File | Size | Format | Contents |
|------|------|--------|----------|
| `cell_feature_matrix.zarr.zip` | 25.1 MB | Zarr (compressed) | Sparse cell × gene count matrix (209,467 cells × 542 features) stored in CSC format |
| `cells.zarr.zip` | 206.2 MB | Zarr (compressed) | Per-cell metadata: centroid coordinates (x, y), cell area, nucleus area, z-level, nucleus count |
| `transcripts.zarr.zip` | 1.18 GB | Zarr (compressed) | Individual transcript records: spatial coordinates (x, y), gene name, and cell assignment ID for every detected RNA molecule |
| `gene_panel.json` | 153 KB | JSON | Panel definition listing all 280+ gene targets with probe names and metadata |

> **Note:** `transcripts.zarr.zip` is large (1.18 GB). It is only loaded in Sections 14 and 15 for spatial MMP11 analysis. All earlier steps use the cell × gene matrix only.

---

## 4. Setup & Installation

### 4.1 Google Drive Setup

Upload all four files to a folder named `xenium_data` in your Google Drive:

```
My Drive/
└── xenium_data/
    ├── cell_feature_matrix.zarr.zip
    ├── cells.zarr.zip
    ├── transcripts.zarr.zip
    └── gene_panel.json
```

### 4.2 Mount Drive in Colab

```python
from google.colab import drive
drive.mount('/content/drive')
```

### 4.3 Install Dependencies

Run the first notebook cell:

```bash
!pip install -q zarr==2.18.3 scanpy anndata shapely geopandas scipy \
             matplotlib seaborn pandas numpy tqdm
!pip install igraph leidenalg
```

> **Why zarr==2.18.3 specifically?** The Xenium zarr files use zarr format v2 conventions. Later versions of zarr (v3+) changed the store API and break compatibility with these files.

### 4.4 Configure Paths

At the top of the notebook, update `DATA_DIR` if your Google Drive folder path differs:

```python
DATA_DIR = '/content/drive/MyDrive/xenium_data'   # ← change if needed
LOCAL    = '/content/xenium_local'                 # temp extraction folder
OUT_DIR  = '/content/xenium_outputs'               # where figures are saved
```

---

## 5. Step-by-Step Code Walkthrough

### Section 0 — Install Dependencies
Installs all required Python packages. Must be run first in every new Colab session.

---

### Section 1 — Imports & Config
Sets up all imports and defines the gene sets used throughout the analysis:

| Variable | Genes | Purpose |
|----------|-------|---------|
| `HK_PAPER` | EEF1G, EEF2, MALAT1, RPLP0 | Paper's selected HK genes for normalization |
| `HK_ONCOTYPE` | ACTB, GAPDH, GUSB, RPLP0, TFRC | Oncotype Dx HK genes (used for comparison) |
| `HK_STABLE` | 15 genes | Full candidate HK gene set plotted in Fig 1A |
| `DEG_GENES` | KIT, TFPI2, PIP, CHI3L1, etc. | Differentially expressed genes (high CV) for Fig 1B |
| `BIOMARKERS` | LDHA, SDC1, SFRP1, PIGR | Tumor grade biomarkers for Figs 1E/1F |

---

### Section 2 — Mount Drive & Extract Zarr Archives
Mounts Google Drive and extracts the three `.zarr.zip` archives to `/content/xenium_local/`. Each archive is only extracted once — if already extracted, it skips to avoid re-doing 1+ GB of decompression.

---

### Section 3 — Load Cell × Gene Matrix → AnnData

**What happens:**  
The cell feature matrix is stored in Zarr as a Compressed Sparse Column (CSC) matrix — efficient for column (gene) access. The code:
1. Reads gene names from `gene_panel.json`
2. Loads `data`, `indices`, `indptr` arrays from the zarr store
3. Constructs a `scipy.sparse.csc_matrix` (genes × cells), then transposes to cells × genes
4. Builds cell IDs from the two-column uint32 `cell_id` array
5. Loads cell summary metadata (centroids, areas, z-level)
6. Assembles everything into an `AnnData` object

**Output:**
```
AnnData object with n_obs=209,467 cells × n_vars=542 genes
obsm: 'spatial'  ← (x, y) centroid coordinates
obs:  cell_centroid_x, cell_centroid_y, cell_area, nucleus_area, ...
```

---

### Section 4 — Load Cell Metadata
Loads and aligns additional cell metadata from `cells.zarr`. Handles potential index mismatches between the two zarr stores by falling back to positional alignment if cell IDs don't match.

---

### Section 5 — QC Filtering

**Paper thresholds applied:**
- Keep cells with **> 10 total transcripts**
- Keep cells with **> 5 unique genes detected**

```
Before QC: 209,467 cells
After QC:  ~175,000–195,000 cells (varies slightly by data version)
```

Also rebuilds the AnnData with the correct feature count (542, not 300), padding unnamed features as `blank_N`. Saves two QC distribution plots:
- Transcript count per cell histogram
- Unique gene count per cell histogram

---

### Section 6 — Normalize → Log1p → PCA → Neighbors → UMAP → Leiden

Standard Scanpy preprocessing pipeline:

| Step | Function | Details |
|------|----------|---------|
| Save raw | `adata.layers['raw']` | Raw counts preserved for CV calculations |
| Normalize | `sc.pp.normalize_total` | Scale each cell to 10,000 total counts |
| Log-transform | `sc.pp.log1p` | Log(x+1) stabilizes variance |
| HVG selection | `sc.pp.highly_variable_genes` | Top 2000 most variable genes |
| PCA | `sc.pp.pca` | 30 components |
| Neighbors | `sc.pp.neighbors` | k-nearest neighbor graph |
| UMAP | `sc.tl.umap` | 2D embedding |
| Clustering | `sc.tl.leiden` | Resolution 0.3 (matches paper) |

---

### Section 7 — Cell Type Annotation
Scores each cell for 7 cell type gene sets using `sc.tl.score_genes`, then assigns each cell the label of the highest-scoring type:

| Cell Type | Marker Genes Used |
|-----------|------------------|
| Tumor | EPCAM, KRT8, KRT18, MKI67, KRT19 |
| Myoepithelial | KRT14, KRT5, ACTA2, MYLK, TP63 |
| T cell | CD3D, CD3E, CD8A, CD4, IL7R |
| B cell | CD79A, MS4A1, CD19 |
| Macrophage | CD68, CSF1R, LYZ |
| Fibroblast | COL1A1, DCN, FAP, VIM |
| Endothelial | PECAM1, VWF, CDH5 |

---

### Section 7b — Assign Tumor Subtypes
The paper uses 4 tissue sections with known subtypes (S2T, S3T, S1B, S2B). Since this notebook uses a single section, tumor cells are sub-clustered using Leiden (resolution 0.5) and the resulting clusters are assigned subtype labels cyclically. This is an approximation — the real biological meaning requires matched multi-section data.

---

### Section 8 — CV Computation Helper

Two utility functions used throughout:

**`get_raw_counts(adata, layer='raw')`**  
Returns the raw count matrix as a dense NumPy array, handling both sparse and dense formats.

**`compute_cv_per_celltype(adata, genes, ...)`**  
For each requested gene and each cell type (with ≥5 cells), computes:
```
CV = std(counts) / mean(counts)
```
Returns a DataFrame of shape `(genes × cell_types)`.

---

### Section 9 — Figure 1A: HK Gene CV Violin Plot

**What it shows:** Distribution of CV values for 15 candidate HK genes, each computed per cell type and shown as a violin + scatter.

**Color coding:**
- 🔵 Teal: Stable HK genes (low CV across cell types) — paper recommendation
- 🟢 Green: Stable in tumor but not universally classified as HK
- 🟠 Orange: Oncotype Dx HK genes

**Interpretation:** Genes like EEF1G, EEF2, MALAT1 should show tight, low violins. Genes like TFRC, GUSB should show wider, higher violins — confirming they are unsuitable as HK references at single-cell resolution.

**Output:** `fig1A_hk_cv_violin.png`

---

### Section 10 — Figure 1B: DEG CV Violin Plot

**What it shows:** CV of 9 differentially expressed genes across tumor subtypes (S2T, S3T, S1B, S2B). These genes are expected to be *highly variable* — that is the whole point. High CV here is a feature, not a bug.

**Color:** All violins in magenta/pink — visually distinguishing them from HK genes.

**Output:** `fig1B_deg_cv_violin.png`

---

### Section 11 — Figure 1C: CV vs Mean Scatter Plot

**What it shows:** Log-log scatter of CV vs. mean expression for all 542 genes in tumor cells. This is a classic noise analysis plot.

**Key insight:** Low-expressing genes naturally have higher CV due to stochastic noise (Poisson statistics). A horizontal and vertical dashed line marks the stable, high-expression zone. Highlighted paper genes (EEF1G, EEF2, MALAT1, RPLP0) should cluster in the low-CV, high-expression quadrant.

**Output:** `fig1C_cv_scatter.png`

---

### Section 12 — Figure 1E: LDHA Raw Expression by Tumor Subtype

**What it shows:** Raw LDHA transcript counts per tumor cell, split by subtype (S2T → S3T → S1B → S2B). Includes jittered scatter, mean line, SEM error bars, and Mann-Whitney U significance brackets.

**Expected trend:** LDHA increases progressively with tumor grade, reflecting the Warburg effect's escalation in more aggressive cancer.

**Output:** `fig1E_LDHA_raw.png`

---

### Section 13 — Figure 1F: LDHA Normalized to HK Genes

**What it shows:** Four panels — LDHA expression normalized to each of EEF1G, RPLP0, GAPDH, and GUSB separately. This is the critical comparison:

- Normalization to **EEF1G** or **RPLP0** should preserve the grade-related LDHA trend
- Normalization to **GAPDH** may flatten or reverse the trend (because GAPDH itself rises with Warburg)
- Normalization to **GUSB** may also distort the trend

**Output:** `fig1F_LDHA_normalized.png`

---

### Section 14 — Figure 2C: Myoepithelial Marker Expression

**What it shows:** Bar chart of KRT14, COL17A1, and LAMC2 expression in myoepithelial cells adjacent to normal ducts vs. tumor (DCIS) ducts.

**Method:** Sub-clusters myoepithelial cells and assigns the first half as "Normal" and the second half as "Tumor-associated" (a rough proxy when working with a single section). Mean log-normalized expression is plotted.

**LAMC2** label appears in red — highlighting it as the key upregulated marker in tumor-associated myoepithelium.

**Output:** `fig2C_myoepithelial_markers.png`

---

### Section 15 — Figures 2F & 2G: MMP11 Periductal Transcript Density

This is the most complex and methodologically novel section.

**Required extra file:** `table.xlsx` (Table S5 from the paper) — contains the X,Y polygon coordinates of 16 manually drawn ROIs from Xenium Explorer (8 around MMP11+ ducts, 8 around MMP11− ducts). Place in `DATA_DIR`.

**Pipeline:**

```
Table S5 ROI polygons
       ↓
Build Shapely Polygon for each ROI
       ↓
Expand outward by 30 µm (PERIPHERY_PX = 30 × 4.706 px/µm)
Subtract inner duct → get periductal ring
       ↓
Load transcripts.zarr (MMP11 + HK genes only)
       ↓
For each ring:
  → Spatial join: find transcripts inside ring
  → Exclude transcripts assigned to duct cells (Tumor, Myoepithelial)
  → Count MMP11 transcripts
  → Count HK gene transcripts (EEF1G, EEF2, MALAT1, RPLP0)
       ↓
Fig 2F: raw_density = MMP11_count / ring_area_µm²
Fig 2G: norm_density = MMP11_count / geomean(HK_counts) / ring_area_µm²
```

**Pixel-to-micron conversion:**  
Xenium uses 0.2125 µm/pixel → 4.706 pixels/µm.  
A 30 µm ring = 141.2 pixels outward expansion.

**Statistical test:** Mann-Whitney U (non-parametric) comparing MMP11+ vs MMP11− ducts.

**Outputs:** `fig2F_MMP11_raw_density.png`, `fig2G_MMP11_normalized_density.png`

---

### Section 16 — Save All Outputs

Saves the final annotated AnnData object and derived tables:

| Output File | Contents |
|-------------|----------|
| `xenium_annotated.h5ad` | Full AnnData with cell type labels, UMAP, Leiden clusters |
| `hk_gene_cv_table.csv` | CV per gene per cell type (for the 15 HK candidate genes) |
| `tableS3_stable_genes.csv` | All genes with minimum CV < 1.5 across any cell type, sorted |

---

## 6. Output Files & Figures

| File | Figure | Description |
|------|--------|-------------|
| `qc_distributions.png` | — | Transcript and gene count distributions post-filter |
| `fig1A_hk_cv_violin.png` | Fig 1A | CV of 15 HK candidate genes across cell types |
| `fig1B_deg_cv_violin.png` | Fig 1B | CV of 9 DEGs across tumor subtypes |
| `fig1C_cv_scatter.png` | Fig 1C | CV vs. mean expression, all genes, tumor cells |
| `fig1E_LDHA_raw.png` | Fig 1E | Raw LDHA Tx/cell across 4 tumor subtypes |
| `fig1F_LDHA_normalized.png` | Fig 1F | LDHA normalized to EEF1G, RPLP0, GAPDH, GUSB |
| `fig2C_myoepithelial_markers.png` | Fig 2C | KRT14, COL17A1, LAMC2 in normal vs tumor-adjacent myoepithelial cells |
| `fig2F_MMP11_raw_density.png` | Fig 2F | MMP11 transcript density (Tx/µm²) in 30 µm periductal ring |
| `fig2G_MMP11_normalized_density.png` | Fig 2G | HK-normalized MMP11 density in periductal ring |
| `xenium_annotated.h5ad` | — | Annotated AnnData object (can be reloaded for further analysis) |
| `hk_gene_cv_table.csv` | — | Gene × cell-type CV table |
| `tableS3_stable_genes.csv` | — | Stable gene list (CV < 1.5) |

---

## 7. Key Tables & Data Summaries

### Cell Count Summary

| Stage | Count |
|-------|-------|
| Raw cells in matrix | 209,467 |
| After QC (>10 Tx, >5 genes) | ~175,000–195,000 |
| Leiden clusters (res=0.3) | ~8–15 |

### Gene Panel Composition

| Category | Count |
|----------|-------|
| Total features in matrix | 542 |
| Named gene targets in panel | 300 |
| Negative controls / blanks | ~242 |
| HK candidate genes | 45 |
| Basement membrane genes | ~30+ |

### HK Gene Selection Summary

| Gene | Category | CV (intra-section) | Recommended |
|------|----------|--------------------|-------------|
| EEF1G | Translation factor | Very low | ✅ Yes |
| EEF2 | Translation factor | Very low | ✅ Yes |
| MALAT1 | lncRNA | Very low | ✅ Yes |
| RPLP0 | Ribosomal protein | Low | ✅ Yes |
| ACTB | Cytoskeletal | Low | ⚠️ Borderline |
| GAPDH | Glycolytic enzyme | High in tumor | ❌ No (Warburg-affected) |
| GUSB | Lysosomal enzyme | Variable | ❌ No |
| TFRC | Transferrin receptor | High variability | ❌ No |

### Pixel-to-Micron Conversion

| Quantity | Value |
|----------|-------|
| Xenium pixel size | 0.2125 µm/pixel |
| Pixels per micron | 4.706 px/µm |
| 30 µm ring in pixels | 141.2 px |
| Ring area (typical) | 50,000–200,000 µm² |

---

## 8. Troubleshooting

### Issue: `zarr.errors.GroupNotFoundError` when loading matrix
**Cause:** The zarr files require zarr v2 API. Newer zarr versions changed the store interface.  
**Fix:** Pin the version: `pip install zarr==2.18.3`

### Issue: Gene count mismatch (300 expected, 542 in matrix)
**Cause:** The gene panel JSON lists 300 named targets, but the matrix has 542 features (including internal controls, blanks, negative probes).  
**Fix:** Already handled in Section 5 — the code pads the name list to 542 with `blank_N` labels.

### Issue: `leidenalg` not found
**Fix:** Run `!pip install leidenalg` then restart the Colab runtime and re-run from Section 1.

### Issue: `No myoepithelial cells found` (Fig 2C skipped)
**Cause:** The annotation in Section 7 depends on which genes are present in the panel. Some panel versions may not include all marker genes.  
**Fix:** Check that `KRT14`, `KRT5`, `ACTA2` are in `adata.var_names`. If fewer than 10 myoepithelial cells are found, the section is automatically skipped with a warning.

### Issue: Fig 2F/2G fails with `FileNotFoundError` for table.xlsx
**Cause:** Table S5 (ROI polygon coordinates) from the paper's supplemental materials must be placed manually in `DATA_DIR`.  
**Fix:** Download Table S5 from the paper's supplemental data and rename it `table.xlsx`.

### Issue: Memory error when loading transcripts.zarr
**Cause:** `transcripts.zarr.zip` (1.18 GB) requires significant RAM.  
**Fix:** The code already filters to only MMP11 and HK transcripts using a boolean mask before loading the full array into memory. Ensure Colab is set to **High RAM** runtime (Runtime → Change runtime type → High-RAM).

### Issue: Tumor subtypes all show as `Unknown`
**Cause:** If the cell annotation step finds zero tumor cells (e.g., none of the tumor markers are in the panel), the subtype assignment defaults to Unknown.  
**Fix:** Verify that `EPCAM`, `KRT8`, `KRT18` are in `adata.var_names`. If not, add alternative epithelial markers.

---

## 9. Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `zarr` | 2.18.3 (pinned) | Read Xenium zarr archives |
| `scanpy` | ≥1.9 | Single-cell analysis pipeline |
| `anndata` | ≥0.9 | AnnData container |
| `numpy` | ≥1.24 | Array operations |
| `pandas` | ≥1.5 | DataFrames |
| `scipy` | ≥1.9 | Sparse matrices, statistics |
| `matplotlib` | ≥3.6 | Plotting |
| `seaborn` | ≥0.12 | Statistical visualization |
| `shapely` | ≥2.0 | Polygon geometry (periductal rings) |
| `geopandas` | ≥0.13 | Spatial join of transcripts to ring polygons |
| `tqdm` | any | Progress bars |
| `igraph` | any | Graph backend for Leiden |
| `leidenalg` | any | Leiden community detection algorithm |
| `openpyxl` | any | Read Table S5 Excel file |

---

## Citation

If you use this notebook or build on this analysis, please cite the original paper:

```bibtex
@article{janesick2025biomarker,
  title   = {Biomarker Quantification in Breast Cancer using Xenium In Situ},
  author  = {Janesick, Amanda and Kravitz, Stephanie N. and Stauffer, Weston 
             and Valencia, Miriam and Taylor, Sarah E.B. and others},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.64898/2025.12.08.692193}
}
```

---

*Notebook author: NUST SINES — Special Topics in Bioinformatics / BI-434*  
*Dataset: 10x Genomics Xenium FFPE Human Breast Biomarkers v4.0.0*
