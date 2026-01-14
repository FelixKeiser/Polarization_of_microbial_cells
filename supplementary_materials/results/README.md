# Results

This directory contains analysis results and outputs generated from the code and data.

## Contents

Typical contents include:
- Statistical test outputs
- Model fitting results
- Summary statistics
- Intermediate analysis results
- Log files from analysis runs

## Organization

```
results/
├── statistical_tests/
│   ├── anova_results.csv
│   └── post_hoc_tests.csv
├── model_outputs/
│   ├── model_parameters.csv
│   └── model_predictions.csv
└── summary_statistics/
    └── descriptive_stats.csv
```

## Documentation

Each results file should include:
- Description of what analysis produced the result
- Date and version of code used
- Parameters used in the analysis
- Link to the code that generated it

## Reproducibility

To ensure reproducibility:
- Document the exact code version used
- Include timestamps
- Note random seeds if applicable
- Record software versions and dependencies

## File Formats

Results should be in standard formats:
- CSV for tabular results
- JSON for structured results
- Plain text for logs and reports
