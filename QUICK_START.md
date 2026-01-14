# Quick Start Guide: Adding Your Supplementary Materials

This repository is now set up to accept supplementary materials for your scientific paper. Follow these simple steps to add your files.

## What's Already Set Up

✅ Directory structure for organizing your materials  
✅ README files in each directory with detailed instructions  
✅ .gitignore configured to handle common file types  
✅ Documentation for citation and usage  

## Quick Steps to Add Your Files

### If you have a .zip file:

1. **Extract your zip file:**
   ```bash
   cd /path/to/this/repository
   unzip your_supplementary_files.zip -d temp_extract
   ```

2. **Organize the files into the correct directories:**
   ```bash
   # Move data files
   mv temp_extract/your_data_files/* supplementary_materials/data/
   
   # Move code files
   mv temp_extract/your_code_files/* supplementary_materials/code/
   
   # Move figures
   mv temp_extract/your_figures/* supplementary_materials/figures/
   
   # Move tables
   mv temp_extract/your_tables/* supplementary_materials/tables/
   
   # Move results
   mv temp_extract/your_results/* supplementary_materials/results/
   ```

3. **Commit and push:**
   ```bash
   git add supplementary_materials/
   git commit -m "Add supplementary materials"
   git push
   ```

### If you want to keep the .zip file as-is:

1. **Place your zip file in the repository:**
   ```bash
   cp your_supplementary_files.zip supplementary_materials.zip
   ```

2. **Update .gitignore to allow this specific file:**
   ```bash
   echo "!supplementary_materials.zip" >> .gitignore
   ```

3. **Commit and push:**
   ```bash
   git add -f supplementary_materials.zip .gitignore
   git commit -m "Add supplementary materials archive"
   git push
   ```

## Directory Structure

```
supplementary_materials/
├── data/           # Raw and processed data files
├── code/           # Analysis scripts and source code
├── figures/        # Supplementary figures
├── tables/         # Supplementary tables
└── results/        # Analysis results and outputs
```

## Important Documents

- **[ADDING_SUPPLEMENTARY_FILES.md](ADDING_SUPPLEMENTARY_FILES.md)** - Detailed guide with multiple methods for adding files
- **[CITATION_TEMPLATE.md](CITATION_TEMPLATE.md)** - Templates for citing these materials in your paper
- **[supplementary_materials/README.md](supplementary_materials/README.md)** - Overview of the supplementary materials structure

## File Size Considerations

- **Small files (< 100 MB)**: Add directly to Git
- **Large files (100 MB - 2 GB)**: Use Git LFS (see ADDING_SUPPLEMENTARY_FILES.md)
- **Very large files (> 2 GB)**: Consider external data repositories (Zenodo, Figshare, etc.)

## Next Steps

1. ✅ Add your supplementary files using one of the methods above
2. ✅ Update README files in each subdirectory to describe your specific files
3. ✅ Add any additional documentation needed
4. ✅ Update CITATION_TEMPLATE.md with your paper's citation information
5. ✅ Test that all your code runs with the uploaded data

## Getting Help

- Check the detailed guides in the repository root
- Review the README files in each subdirectory
- See [ADDING_SUPPLEMENTARY_FILES.md](ADDING_SUPPLEMENTARY_FILES.md) for common issues and solutions

## Checklist Before Submission

- [ ] All data files are in `supplementary_materials/data/`
- [ ] All code files are in `supplementary_materials/code/`
- [ ] All figures are in `supplementary_materials/figures/`
- [ ] All tables are in `supplementary_materials/tables/`
- [ ] Each subdirectory has an updated README describing the contents
- [ ] Code has been tested and runs correctly
- [ ] File naming is consistent and descriptive
- [ ] Large files are handled appropriately (Git LFS or external storage)
- [ ] Citation information has been updated
- [ ] All changes are committed and pushed to GitHub

---

**Ready to start?** Begin by extracting your .zip file and organizing the contents into the appropriate directories!
