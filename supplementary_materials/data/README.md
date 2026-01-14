# Data

This directory contains data files for the study on polarization of microbial cells.

## Subdirectories

- `raw/` - Original, unmodified data files as collected from experiments
- `processed/` - Cleaned and processed data ready for analysis
- `metadata/` - Data dictionaries, column descriptions, and other metadata

## File Naming Convention

Please use descriptive names following this pattern:
- `YYYY-MM-DD_experiment_description.csv`
- `dataset_name_processed.csv`

## Data Format

Data files should preferably be in open formats:
- CSV for tabular data
- JSON for structured data
- Plain text for simple data

## Documentation

Each data file should be accompanied by documentation describing:
- Variables/columns and their meanings
- Units of measurement
- Experimental conditions
- Date of collection
- Any processing steps applied

## Example Structure

```
data/
├── raw/
│   ├── 2024-01-15_cell_tracking_exp1.csv
│   └── 2024-01-20_cell_tracking_exp2.csv
├── processed/
│   ├── combined_cell_tracking.csv
│   └── normalized_cell_data.csv
└── metadata/
    └── data_dictionary.csv
```
