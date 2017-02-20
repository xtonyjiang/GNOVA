# Under Construction

## Requirements
1. Python 2.7
2. numpy
3. scipy
4. pandas
5. sklearn
6. bitarray

## Example Usage

```
python pipeline.py data/CD.sumstats.gz data/UC.sumstats.gz --N1 27726 --N2 28738 --bfile data/bim/eur_chr@_SNPmaf5 --annot data/annot/func.@.txt --out results.txt
```

## Credits
LD score calculation adapted from LDSC. See
[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)
