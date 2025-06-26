# D-Sites: A High-Throughput Computational Tool for Predicting Protein Binding Sites in Bacterial Genomes

![Python](https://img.shields.io/badge/Python-3.12-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-Linux--Windows--Unix-lightgrey)

---

## Overview

**D-Sites** is a modular, high-throughput computational pipeline for the **genome-wide prediction of DNA-binding protein interaction sites** in microbial genomes. It combines motif matching, DNA structural feature profiling, and palindromic scoring into a unified scoring framework to identify potential protein–DNA binding regions with biological relevance.

>  **Live Web Application**:  
> Try the tool online at:  
> 🔗 **[https://sinr-backend.onrender.com](https://sinr-backend.onrender.com)**

---

##  Key Features

- Predicts binding sites for any DNA-binding protein with known motif(s)
- Incorporates DNA shape descriptors: Minor Groove Width, Propeller Twist, and Roll
- ♻Efficient parallel processing of complete microbial genomes
- Outputs interactive plots, scores, and binding coordinates
- Lightweight, no GPU required
- Tested on *Salmonella enterica* using SinR as a case study

---

## Repository Contents

- `protein_binding_site_prediction.py` – Core pipeline
- `requirements.txt` – Python dependencies
- `example_data/` – Sample genome and protein FASTA files
- `output/` – Example CSV and plot outputs

---

## Installation

### Requirements

- Python 3.10+
- Modules: `numpy`, `pandas`, `biopython`, `matplotlib`, `seaborn`, `tqdm`

```bash
pip install -r requirements.txt
```
