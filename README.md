# GNOVA

GNOVA (GeNetic cOVariance Analyzer), a principled framework to estimate annotation-stratified genetic covariance using GWAS summary statistics.

## Requirements
1. Python 2.7
2. numpy
3. scipy
4. pandas
5. sklearn
6. bitarray

## Tutorial

Suppose you would like to calculate genetic covariance between Crohn's Disease and Ulcerative Colitis. We'll need a few types of files:

- **Summary statistics files:** You can get your own GWAS summary statistics files for these two diseases [here](https://www.ibdgenetics.org). We assume that the files are in the standard format that `ldsc` understands. If not, make sure to run them through the included `munge_sumstats.py` file, or use the one included in `ldsc` (see [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics) for instructions).

- **Plink bfiles:** These are files .bed/.bim/.fam format. You can download some that we have prepared for you [here](http://genocanyon.med.yale.edu/GNOVAFiles/genotype_1KG_eur_SNPmaf5.tar.gz). These files are from the 1000 Genomes Project, with rare variants (MAF < 5%) filtered out.

- **Anotation files:** These are only necessary if you are doing annotation-stratified analysis. You can download some example annotation files that we have prepared for you [here](http://genocanyon.med.yale.edu/GNOVAFiles/annotations.tar.gz). There are three different types of annotation files here:
    - **annot_GenoCanyon:** Predicted functional and non-functional SNPs based on GenoCanyon annotations.
    - **annot_MAF:** Minor allele frequency quartiles based on samples from the 1000 Genomes Project with European ancestry.
    - **annot_GSplus7tissues:** Seven broadly-defined tissue categories based on GenoSkyline-Plus annotations.

More details about these supplied files can be found in in the [GNOVA manuscript](http://biorxiv.org/content/early/2017/03/07/114561).

We may run the following command:

```
python gnova.py data/CD.sumstats.gz data/UC.sumstats.gz \
--N1 27726 \
--N2 28738 \
--bfile data/bfiles/eur_chr@_SNPmaf5 \
--annot data/annot/func.@.txt \
--out results.txt
```
### Explanation of Command-Line Arguments

- The first two arguments, `data/CD.sumstats.gz` and `data/UC.sumstats.gz`, denote the locations of the first and second summary statistics files. These files may be compressed using gzip, bz2, zip, xz, or not compressed at all. The program will infer the compression method if the files end with .gz, .bz2, .zip, xz, respectively. As previously mentioned, we assume that the files are in the standard format that `ldsc` understands.

- The `N1` and `N2` arguments (optional) denote the sample sizes of the summary statistics files. If they are not provided, they will be inferred from the summary statistics files.

- The `bfile` argument denotes the prefix of the `.bed/.bim/.fam` genotypic data files. Note the '@', which denotes a wildcard character that GNOVA will be replace with 1-22. Alternatively, if you would only like analysis on one chromosome, you can just specify that one bfile.

- The `annot` argument (optional) denotes the location of the annotation files if doing annotation-stratified analysis. **We assume that for each chromsome that we are doing analysis on, there is a corresponding whitespace-delimited annotation file for that chromosome, such that if there are n rows in the bim file for chromosome 1, there are n+1 rows for the corresponding annotation file (the annotation file should have an extra row denoting the names of the annotations).**

- The `out` flag denotes the file location for the results to be outputted to.

### Additional Command-Line Arguments

Here is an explanation of the other command-line arguments that weren't shown in the example:

- `--save-ld`: Running GNOVA with `--save-ld ldscores` will save intermediate LD score calculations to `ldscores.csv.gz`.

- `--use-ld`: Running GNOVA with `--use-ld ldscores` will skip the LD score calculation step and instead load the intermediate calculations from `ldscore.csv.gz` that were generated with `--save-ld`.

### Explanation of Output
The output will be a whitespace-delimited text file, with the rows corresponding to different annotations and the columns as such:

- `rho:` The genetic covariance estimate.
- `rho_corrected:` The genetic covariance estimate with sample overlap correction.
- `pvalue:` The p-value from the statistical test for genetic covariance.
- `pvalue_corrected:` The p-value from the statistical test for genetic covariance with sample overlap correction.
- `corr`: The genetic correlation estimate.
- `corr_corrected`: The genetic correlation estimate with sample overlap correction.
- `h2_1`: The heritability estimate for the first trait.
- `h2_2`: The heritability estimate for the second trait.

NOTE: When functional annotations are present, the true heritability in each annotation category may be small. Although methods for estimating annotation-stratified heritability exist, they may provide unstable, in many cases negative heritability estimates, especially when a number of annotation categories are related to the repressed or non-functional genome. GNOVA ignores negative hertiability estimates, leaving the correlation estimates as 'NaN'. So, we recommend the users to focus on genetic covariance instead of genetic correlation when performing annotation-stratified analysis.

## Credits
Those using the GNOVA software should cite:

[Lu, et al. A powerful approach to estimating annotation-stratified genetic covariance using GWAS summary statistics. bioRxiv, 2016.](http://biorxiv.org/content/early/2017/03/07/114561)

The LD score calculation is adapted from `ldsc`. See
[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)
