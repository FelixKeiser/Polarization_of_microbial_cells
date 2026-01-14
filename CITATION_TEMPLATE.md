# Supplementary Materials Citation Template

This document provides templates for referencing the supplementary materials in your manuscript.

## In-Text Citations

### For specific supplementary figures:
```
...as shown in Supplementary Figure S1 (see Supplementary Materials).
```

### For specific supplementary tables:
```
...detailed results are provided in Supplementary Table S1 (see Supplementary Materials).
```

### For code:
```
...analysis was performed using custom Python scripts available in the Supplementary Materials (supplementary_materials/code/).
```

### For data:
```
...raw data are available in the Supplementary Materials (supplementary_materials/data/).
```

## Methods Section

### Data Availability Statement
```
Data Availability

Raw and processed data supporting the findings of this study are available in 
the supplementary materials repository at [GitHub repository URL]. The repository 
includes:
- Raw experimental data (supplementary_materials/data/raw/)
- Processed datasets (supplementary_materials/data/processed/)
- Analysis code (supplementary_materials/code/)
- Supplementary figures (supplementary_materials/figures/)
- Supplementary tables (supplementary_materials/tables/)

[If using external storage:]
Large data files (>100MB) are available at [DOI/URL].
```

### Code Availability Statement
```
Code Availability

All code used for data analysis and figure generation is available in the 
supplementary materials repository at [GitHub repository URL] under the 
[LICENSE TYPE] license. The code includes:
- Data preprocessing scripts (supplementary_materials/code/preprocessing/)
- Statistical analysis scripts (supplementary_materials/code/analysis/)
- Visualization code (supplementary_materials/code/visualization/)

Requirements and usage instructions are provided in the repository README.
```

## Supplementary Materials Section (for manuscript)

```
Supplementary Materials

The following supplementary materials are available at [GitHub repository URL]:

Supplementary Figures:
- Figure S1: [Brief description]
- Figure S2: [Brief description]

Supplementary Tables:
- Table S1: [Brief description]
- Table S2: [Brief description]

Supplementary Data:
- Dataset S1: [Brief description]
- Dataset S2: [Brief description]

Supplementary Code:
- Script S1: [Brief description and filename]
- Script S2: [Brief description and filename]
```

## README Section (to add to repository README)

```markdown
## Citation

If you use these supplementary materials, please cite:

[Author names]. ([Year]). [Paper title]. [Journal name], [Volume]([Issue]), 
[Page range]. DOI: [DOI]

BibTeX:
```bibtex
@article{authorYYYY,
  title={Paper title},
  author={Author, First and Author, Second},
  journal={Journal name},
  volume={XX},
  number={X},
  pages={XXX--XXX},
  year={YYYY},
  publisher={Publisher},
  doi={10.xxxx/xxxxx}
}
```

## For Grant Proposals / Publications

When describing data management plans:

```
All research data and analysis code will be made publicly available via a 
GitHub repository [URL]. The repository will include:
1. Raw and processed data in open formats (CSV, JSON)
2. Analysis scripts with documentation
3. Supplementary figures in publication-ready formats
4. Detailed README files describing the repository structure and usage

Data will be preserved for at least [X] years following publication and will 
be licensed under [LICENSE NAME] to facilitate reuse and reproducibility.
```

## Journal-Specific Requirements

### For journals requiring structured supplementary materials:

Create a file `SUPPLEMENTARY_FILE_LIST.md`:

```markdown
# Supplementary File List

## Supplementary Figures
1. **supplementary_fig_S1.pdf** - Cell polarization under different conditions
2. **supplementary_fig_S2.pdf** - Time-lapse analysis results

## Supplementary Tables
1. **supplementary_table_S1.csv** - Complete statistical analysis results
2. **supplementary_table_S2.csv** - Cell measurement data

## Supplementary Data
1. **raw_cell_tracking.csv** - Raw cell tracking data
2. **processed_measurements.csv** - Processed cell measurements

## Supplementary Code
1. **analysis_pipeline.py** - Main analysis pipeline
2. **visualization.R** - Figure generation scripts
```

## DOI Assignment

For long-term preservation, consider obtaining a DOI through:
- **Zenodo**: Automatically assigns DOIs to GitHub releases
- **figshare**: Provides DOIs for datasets and code
- **Dryad**: Specialized for scientific data

To create a Zenodo DOI from this GitHub repository:
1. Go to https://zenodo.org/
2. Log in with your GitHub account
3. Enable the repository in your Zenodo GitHub settings
4. Create a release in GitHub
5. Zenodo will automatically create a DOI

Then update your citations with the DOI:
```
DOI: 10.5281/zenodo.XXXXXXX
```
