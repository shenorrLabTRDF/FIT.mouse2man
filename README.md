# FIT.mouse2man

Cross-species differences form barriers to translational research that ultimately hinder the success of clinical trials. 
Yet systematic incorporation of the wealth of knowledge on species differences in the interpretation of animal model data 
is lacking. Here we present FIT (Found In Translation) a data-driven statistical methodology which leverages public domain 
gene expression data to predict from the results of a mouse experiment genes expected to be altered in the equivalent human phenotype.

## Getting Started

First you will need to download the package from github.
Start by installing the devtools package from CRAN and load it:
```
install.packages("devtools")
library(devtools)
```

Install FIT_mouse2man and load it:
```
require(devtools)
devtools::install_github('shenorrLab/FIT.mouse2man')
```

### Checking applicability of the model to your data
FIT is able to improve the similarity to human data in a wide variety of diseases, yet not in all of them. FIT's performance relies, among other features, on the mouse input data. To allow users to know in advance whether FIT is likely or unlikely to improve the similarity of their mouse experiment to human data, we created an SVM classifier that predicts FIT performance. The classifier is accurate on avrage in 80% of tested cases.
The classifier is given:
 - The mouse data file in CSV format (see below).
 - q-value threshold - the q-value by which differential genes are going to be defined (Default: 0.1).
 - Fold-change threshold - the fold-change by which differential genes are going to be defined. The fold-change is defined as fraction of genes from the top. For example, fold-chnage=0.15 denotes 15% of the genes with highest absolute fold-chnage (Default: 0.15).
 The latter two parameters are not used in any analysis done by the package, but affect the classifier response space.
 
Example: 
```
> RunClassifier("MyMouseData.csv")
The mouse data contains 11545 genes.
The classifier can be based on 4957 genes.
The current run will be based on 4957 genes (intersection between the current data and the classifier set of genes.

*************************  Classifier prediction  **************************
It is unlikely FIT will be able to improve this dataset.
****************************************************************************

See the performance results of the classifier to identify the performance of the classifier in the selected set of theesholds (Fold-change=0.15, q-value=0.1)
```

### Input data
The input data for the model is a mouse gene expression dataset in CSV (comma separated values) format.

#For microarray data:
The rows of the table represent genes and the columns are samples. The values should be log-transformed expression values. 
The first column should contain the gene names, as Entrez IDs. The first row should contain the sample names.
The dataset should contain at least 3 control samples and at least 3 disease samples. 
The sample names should include an annotation of which sample is control and which is disease, in the following way: disease sample names should start with "d_" and control sample names should start with "c_". 
If one of the following criteria regarding the file format is not met, an error message will be shown: 
-	Sample names start with either "c_" or "d_". 
-	There are at least 3 control and least 3 disease samples. 
-	The gene names are mouse Entrez IDs. 

#For RNAseq data:
The CSV file should contain 2 columns: the first is mosue Entrez gene IDs, and the second is effect size per gene. The effect size can be aquired by using a RNAseq data analysis software, for exmaple Kallisto and Sleuth (https://pachterlab.github.io/kallisto/, https://pachterlab.github.io/sleuth/).

Sample files can be obtained by the function GetSampleData().

## Example
To run the FIT pipeline, use the FIT() function:
```
res = FIT("MyMouseData.csv", "microarray") # Runs the FIT pipeline and output prediction per gene
```

The output is sorted by the absolute value of FIT prediction:
```
res[1:3,]
      Mouse.Entrez Human.Entrez Mouse.symbol Human.symbol Mouse_FoldChange Mouse_Ztest FIT_prediction FIT_percentile UpDown   Conf_low
2283         14265         2332         Fmr1         FMR1      -2.43926817  -25.299665     -14.624135         100.00      - -18.446825
11355        67213        54918        Cmtm6        CMTM6       0.16115661    4.993052       5.564122          99.99      +   4.532928
10006        56401        64175       Lepre1       LEPRE1       0.10606464    4.168578       4.236725          99.99      +   2.657927
       Conf_high   CI_size CI_percentile
2283  -10.106493 8.3403321        100.00
11355   5.993052 1.4601242         96.65
10006   5.168578 2.5106516         99.54
```

The results allow to focus on genes that are relevant to the human disease. The recommendation is to focus on genes with a high absolute predicted value and low confidence interval sizes.

## Authors

* **Rachelly Normand, Wenfei Du** 
