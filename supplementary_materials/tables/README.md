# Tables

This directory contains supplementary tables referenced in the manuscript.

## File Formats

Tables should be provided in accessible formats:
- **CSV** - Comma-separated values (preferred for data)
- **Excel (.xlsx)** - For complex tables with formatting
- **Markdown (.md)** - For simple tables included in documentation

## Naming Convention

- `supplementary_table_S1.csv`
- `supplementary_table_S2.xlsx`

## Table Documentation

Each table file should include:
- A header row with clear column names
- A separate README or comments explaining:
  - What the table contains
  - Column definitions
  - Units of measurement
  - Any abbreviations used

## Example Structure

```
tables/
├── supplementary_table_S1.csv
├── supplementary_table_S1_description.md
├── supplementary_table_S2.csv
└── supplementary_table_S2_description.md
```

## Best Practices

- Use consistent column naming across all tables
- Include units in column headers (e.g., "Length (μm)")
- Use standard missing data indicators (e.g., NA, NULL)
- Avoid special characters in column names
- Save in UTF-8 encoding
