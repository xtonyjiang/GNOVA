# GNOVA

GNOVA (GeNetic cOVariance Analyzer), a principled framework to estimate annotation-stratified genetic covariance using GWAS summary statistics.

## Requirements
1. Python 2.7
2. numpy
3. scipy
4. pandas
5. sklearn
6. bitarray

## Example Usage

Suppose you would like to calculate genetic covariance between Crohn's Disease and Ulcerative Colitis. We may run the following command:

```
python gnova.py data/CD.sumstats.gz data/UC.sumstats.gz \
--N1 27726 \
--N2 28738 \
--bfile data/bim/eur_chr@_SNPmaf5 \
--annot data/annot/func.@.txt \
--out results.txt
```
### Explanation of Command-Line Arguments

- The first two arguments, `data/CD.sumstats.gz` and `data/UC.sumstats.gz`, denote the locations of the first and second summary statistics files. These files may be compressed using gzip, bz2, zip, xz, or not compressed at all. The program will infer the compression method if the files end with .gz, .bz2, .zip, xz, respectively. We assume that the files are in the standard format that `ldsc` understands. If not, make sure to run them through the included `munge_sumstats.py` file, or use the one included in `ldsc` (see [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics) for instructions)

- The `N1` and `N2` arguments (optional) denote the sample sizes of the summary statistics files. If they are not provided, they will be inferred from the summary statistics files.

- The `bfile` argument denotes the prefix of the `.bed/.bim/.fam` genotypic data files. The '@' denotes a wildcard character that will be replaced by 1-22, for multi-chromosome analysis.

- The `annot` argument (optional) denotes the location of the annotation files if doing annotation-stratified analysis. **We assume that for each chromsome that we are doing analysis on, there is a corresponding whitespace-delimited annotation file for that chromosome, such that if there are n rows in the bim file for chromosome 1, there are n+1 rows for the corresponding annotation file (the annotation file should have an extra row denoting the names of the annotations)**

- The `out` flag denotes the file location for the results to be outputted to.

## Credits
LD score calculation adapted from LDSC. See
[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)
