# FIT.mouse2man

Cross-species differences form barriers to translational research that ultimately hinder the success of clinical trials. 
Yet systematic incorporation of the wealth of knowledge on species differences in the interpretation of animal model data 
is lacking. Here we present FIT (Found In Translation) a data-driven statistical methodology which leverages public domain 
gene expression data to predict from the results of a mouse experiment genes expected to be altered in the equivalent human phenotype.

## Getting Started

First you will need to download the package from Github.
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
FIT is able to improve the similarity to human data in a wide variety of diseases, yet not in all of them. FIT's performance relies, among other features, on the mouse input data. To allow users to know in advance whether FIT is likely or unlikely to improve the similarity of their mouse experiment to human data, we created an SVM classifier that predicts FIT performance. The classifier is accurate on average in 80% of tested cases.
The classifier is given:
 - The mouse data file in CSV format (see below).
 - q-value threshold - the q-value by which differential genes are going to be defined (Default: 0.1).
 - Fold-change threshold - the fold-change by which differential genes are going to be defined. The fold-change is defined as fraction of genes from the top. For example, fold-change=0.15 denotes 15% of the genes with highest absolute fold-change (Default: 0.15).
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

See the performance results of the classifier to identify the performance of the classifier in the selected set of thresholds (Fold-change=0.15, q-value=0.1)
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
The CSV file should contain 2 columns: the first is mouse Entrez gene IDs, and the second is effect size per gene. The effect size can be acquired by using a RNAseq data analysis software, for example Kallisto and Sleuth (https://pachterlab.github.io/kallisto/, https://pachterlab.github.io/sleuth/).

Sample files can be obtained by the following commands:
```
data(microarray_sample)
data(RNAseq_sample)
```

### Example
To run the FIT pipeline, use the FIT() function:
```
res = FIT("MyMouseData.csv", "microarray") # Runs the FIT pipeline which computes prediction of the human effect-size per gene
```

The output is sorted by the absolute value of FIT prediction:
```
res[1:3,]
      Mouse.Entrez Human.Entrez Mouse.symbol Human.symbol Description Mouse_FoldChange Mouse_EffectSize FIT_prediction FIT_percentile UpDown   CI_low  CI_high  CI_size CI_percentile
7878        56277        55076      Tmem45a      TMEM45A     Tmem45a        1.1519944         3.574920       4.403494         100.00      + 3.347403 4.574920 1.227516         96.41
5348        22004         7169         Tpm2         TPM2        Tpm2        0.3530030         5.086331       4.172559          99.99      + 2.890907 5.274678 2.383770         99.83
3911        19128         5627        Pros1        PROS1       Pros1        0.7826637         3.813684       4.129875          99.98      + 3.053448 4.813684 1.760236         99.14
```

The results allow to focus on genes that are relevant to the human disease. The recommendation is to focus on genes with a high absolute predicted value and low confidence interval sizes.


### Training data
The training data which is used by the model is available by:
```
data(AllData_V2.0)
```



## Authors

* **Rachelly Normand, Wenfei Du** 
