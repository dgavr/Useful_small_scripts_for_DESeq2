# Useful small python scripts for dealing with DESeq2 differential expression 


## Annotate differential expression or any other table with another list 

Useful if using Ensembl gene models and would prefer to have gene names. If using Ensembl biomart to fetch geneid, "Associated gene name" and "Gene Description" in that order, allows to annotated in the first three columns these fields, making it easier to make sense of differential expression datasets. 

```{r eval = FALSE}
annotate_count_table2.py zv_mart_export.txt diff_exp_file.txt
```
 
This will generate a new file called annotated_diff_exp_file.txt 

## Trying different combinations of DESeq2  

If faced with non-reproducible RNA-seq datasets and do not know which ones to remove and do not want to try all options by hand - script that runs different combinations and writes then to file. Results can then be screened for positive control genes by simply using grep in the command line. 

```{r eval = FALSE}
combinatorials2.py [Experimental conditions separated by comma] [Names of samples of treated condition separated by comma] [Names of control samples separated by comma] [Name of count table file]
```
