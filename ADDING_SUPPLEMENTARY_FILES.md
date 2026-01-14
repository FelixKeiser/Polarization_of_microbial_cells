# How to Add Your Supplementary Files

This guide explains how to add your supplementary materials from a .zip file to this repository.

## Method 1: Extract and Add to Git (Recommended for Small Files)

If your supplementary files are relatively small (< 100 MB total):

1. Extract your .zip file:
   ```bash
   unzip your_supplementary_files.zip -d supplementary_materials/
   ```

2. Review the extracted files and organize them into the appropriate directories:
   - Data files → `supplementary_materials/data/`
   - Code files → `supplementary_materials/code/`
   - Figures → `supplementary_materials/figures/`
   - Tables → `supplementary_materials/tables/`
   - Results → `supplementary_materials/results/`

3. Add and commit the files:
   ```bash
   git add supplementary_materials/
   git commit -m "Add supplementary materials"
   git push
   ```

## Method 2: Using Git LFS (For Large Files)

If your files are large (> 100 MB), consider using Git Large File Storage:

1. Install Git LFS:
   ```bash
   git lfs install
   ```

2. Track large file types:
   ```bash
   git lfs track "*.zip"
   git lfs track "*.hdf5"
   git lfs track "*.mat"
   # Add other large file extensions as needed
   ```

3. Add and commit your files:
   ```bash
   git add .gitattributes
   git add supplementary_materials/
   git commit -m "Add large supplementary materials using LFS"
   git push
   ```

## Method 3: Keep the .zip File

If you prefer to keep the supplementary materials as a .zip file:

1. Place the .zip file in the repository root:
   ```bash
   cp your_supplementary_files.zip supplementary_materials.zip
   ```

2. Update `.gitignore` to allow this specific .zip file:
   ```bash
   # In .gitignore, add:
   !supplementary_materials.zip
   ```

3. Add and commit:
   ```bash
   git add -f supplementary_materials.zip
   git commit -m "Add supplementary materials archive"
   git push
   ```

## Method 4: External Storage (For Very Large Files)

For very large datasets (> 1 GB):

1. Upload to a data repository:
   - Zenodo (https://zenodo.org/)
   - Figshare (https://figshare.com/)
   - Dryad (https://datadryad.org/)
   - OSF (https://osf.io/)

2. Add a link in the README:
   ```markdown
   ## Large Data Files
   
   Large supplementary data files are available at:
   - [Zenodo DOI: 10.5281/zenodo.XXXXXX](https://zenodo.org/record/XXXXXX)
   ```

## Best Practices

- **Document everything**: Add README files to each subdirectory
- **Use consistent naming**: Follow a clear naming convention
- **Include metadata**: Add data dictionaries and descriptions
- **Version control**: Commit changes incrementally
- **Check file sizes**: Be mindful of repository size limits
- **Test your code**: Ensure scripts run on the uploaded data
- **Add licenses**: Specify usage rights for your materials

## File Organization Tips

```
supplementary_materials/
├── data/
│   ├── raw/                  # Original, unprocessed data
│   ├── processed/            # Cleaned and processed data
│   └── metadata/             # Data descriptions
├── code/
│   ├── analysis/             # Analysis scripts
│   ├── preprocessing/        # Data cleaning scripts
│   └── visualization/        # Plotting scripts
├── figures/
│   ├── supplementary_fig_1/  # One directory per figure
│   └── supplementary_fig_2/
├── tables/
│   └── supplementary_table_1.csv
└── results/
    ├── statistical_tests/
    └── model_outputs/
```

## Need Help?

- Check GitHub's documentation on [file size limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files)
- Learn more about [Git LFS](https://git-lfs.github.com/)
- Review [data repository options](https://www.nature.com/sdata/policies/repositories)
