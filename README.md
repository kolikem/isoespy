# Isoespy

**Isoespy** is a visualization and functional annotation tool for transcript isoforms identified from long-read transcriptome data. It assists in the interpretation of isoform-level transcriptomic models, which are identified and quantified from RNA-seq data, by supporting differential expression analysis, novel isoform ORF prediction, and integration of external functional predictions. The tool provides intuitive visualizations of gene expression states, based on these processed data.

---

## Features

- Differential expression visualization between two sample groups (e.g., tumor vs. normal)
- Visualization of:
  - CDS (coding sequence) predictions
  - Pfam domains
  - Signal peptides
  - NMD-susceptible isoforms
- Visualizations are based on transcript-level GTF files

---

## Installation

Isoespy can be installed via PyPI:

```bash
pip install isoespy
```

or

```bash
python -m pip install isoespy
```

*Requires Python 3.8+*

---

## Commands Overview

Isoespy provides six CLI commands:

### 1. `isoespy-orfpred`

Predicts ORFs for unannotated transcripts using CPAT and TransDecoder.

```bash
isoespy-orfpred --workdir <dir> --hexamer <Path to hexamer table> --model <Path to model> \
                --fasta <nt fasta> --cpatoutprefix <prefix> --transdecoder_path <Path to TransDecoder> \
                --gtf <GTF file> --meta <metadata>
```

---

### 2. `isoespy-makefa`

Extracts transcript nucleotide or amino acid sequences from GTF and genome reference.

```bash
isoespy-makefa --gtf <GTF file> --genome <genome reference> --meta <metadata> \
               --feature <exon/CDS> --type <nucleotide/amino_acid>
```

---

### 3. `isoespy-edger`

Performs DE analysis using edgeR with a prepared R script and sample metadata.

```bash
isoespy-edger --r_script <R script> --meta_data <metadata> --count_data <expression count> \
              --output_de <DE output> --output_cpm <CPM output>
```

---

### 4. `isoespy-makegtf`

Merges external annotation TSV into a GTF file for visualization.

```bash
isoespy-makegtf --gtfprep_file <input TSV> --gtf_file <original GTF> --meta_file <metadata> \
                --output <output GTF> --frame_file <frame info file from makefa>
```

---

### 5. `isoespy-de`

Visualizes DE results for a gene with intron compression and outlier control.

```bash
isoespy-de --gene_name <gene> --gtf_data <GTF> --expression_data <expression> \
           --det_data <DE result> --meta_data <metadata> --compress_introns <int size> \
           --show_outliers <True/False>
```

---

### 6. `isoespy-ff`

Visualizes functional annotations (domains, signals, etc.) for isoforms.

```bash
isoespy-ff --gene_name <gene> --gtf_data <GTF> --meta_data <metadata> \
           --compress_introns <int size> --annotation <annotation TSV>
```

---

## Requirements

### Files

- Reference genome FASTA and index
- CPAT pretrained hexamer and model files

### External Tools

- [CPAT](https://github.com/huangyh09/CPAT)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)

### Python Libraries

- `pysam`
- `matplotlib`
- `pandas`
- `seaborn`
- `numpy`

### R Libraries

- `edgeR`

Users are responsible for preparing the required files and ensuring `edgeR` is installed in their R environment.

---

## License

MIT License

Developed by **Ko Ikemoto**
