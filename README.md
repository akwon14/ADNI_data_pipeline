# ADNI Data Pipeline

Data wrangling pipeline for ADNI amyloid and tau PET CSV files.

## Overview
This pipeline:
- Merges amyloid and tau PET data with ADNI registry visits
- Applies PET QC filtering
- Matches scans by nearest visit date (by RID)
- Integrates ADNIMERGE clinical variables
- Maps APOE genotype to APOE4 status
- Aligns FreeSurfer cortical SUVR regions

## Requirements
- Python >= 3.9
- pandas
- numpy

## Usage
1. Download required ADNI CSV files (not included)
2. Update file paths in `adni_pet_pipeline.py`
3. Run:
   ```bash
   python adni_pet_pipeline.py
