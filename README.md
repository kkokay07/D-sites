# D-Sites: A Computationally efficient tool for predicting protein binding sites in microbial genomes

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
## About the Author

**Dr. Kanaka K. K., PhD, ARS**  
Scientist  
School of Bioinformatics  
[ICAR-Indian Institute of Agricultural Biotechnology, Ranchi](https://iiab.icar.gov.in/)
> [Be like IIAB!:](https://www.researchgate.net/publication/379512649_ICAR-IIAB_Annual_Report-_2023) IIAB is like yogic center where all the sciences (Plant, Animal, Aquatic,Mibrobiology, IT) meet to address emerging issues in food production.

## Spy on me
- [Google Scholar](https://scholar.google.com/citations?hl=en&user=0dQ7Sf8AAAAJ&view_op=list_works&sortby=pubdate)
- [GitHub: kkokay07](https://github.com/kkokay07)
- [ResearchGate](https://www.researchgate.net/profile/Kanaka-K-K/research)
- [Institute Website](https://iiab.icar.gov.in/staff/dr-kanaka-k-k/)
