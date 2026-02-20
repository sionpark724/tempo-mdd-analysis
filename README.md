# Circadian Transcriptomic Analysis of Human MDD Brain Tissue

## Overview
This repository contains the Python scripts and data analysis pipeline for investigating chronobiological disruptions in Major Depressive Disorder (MDD). The analysis applies the **TEMPO** algorithm (Auerbach et al., 2022) to estimate circadian phases and amplitudes from human postmortem bulk RNA-seq data.

## Dataset
* **Source:** GEO Accession `GSE102556` (Labonte et al., 2017)
* **Samples:** 261 human postmortem brain tissue samples.
* **Brain Regions Analyzed:** OFC, dlPFC, vmPFC, aINS, NAc, vSUB.

## Methodology
* **Algorithm:** TEMPO (Python-based)
* **Reference Gene:** *ARNTL* was set to 0h for temporal alignment.
* **Core Clock Genes Prior:** Custom basic prior list (`core_clock_and_ubiq_acrophase_prior.csv`).
* **Random Seed:** Fixed to `2026` for reproducibility.

## Key Findings
1. **Global Rhythm Dampening:** Core clock genes exhibited consistent rhythm dampening (reduced amplitude) in the MDD cohort across multiple brain regions.
2. **Universal Disrupted Hub Genes:** Identified ***NFIL3*, *CIART*, and *CSNK1A1*** as the most severely altered genes, showing significant phase shifts and amplitude vanishing across all six regions.

## Next Steps
* Overcome tissue-specific sample size limitations by pooling data to conduct robust sex-specific circadian analyses.
* Target specific stress susceptibility hub genes (e.g., *DUSP6*) for rhythmic alteration analysis.

## Requirements
* Python 3.x
* tempo (https://github.com/bauerbach95/tempo)
* pandas, numpy, matplotlib, seaborn
