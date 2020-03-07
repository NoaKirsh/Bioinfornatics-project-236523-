#  👨‍🔬 Bioinfornatics-project - 236523 👩‍🔬

Studying the influence of sequencing coverage / replicates on differential analysis.

## ❓ Main questions:

In our project we will attempt to answer the following questions:
*	Primary: What is the ideal balance between coverage and replicates for a given budget?

*	Ideal balance depending on amount of budget

* Replicates Origin: Tradeoff between number of individuals

and number of replicates from each one.

## 💼 Workflow:

* Load the data.

* Removing all the low quality data based on the correlation to the ERCC genes.

* Choosing random samples for a given number of mice and percentage of samples.

* Coverage - Generating random gene samplings from the data, based on their probability to be chosen in the “real world” for different.

* Creating rawcounts and metadata files for deseq, including all the genes and samples for a  given number of mice and coverage percentage.

* Repeating the process for many combinations of different samples and different coverage percentages in order to simulate a variety of budgets (total number of reads) and replicate origin combinations.

* DESEQ - running deseq to analyse sequencing data and test for differential Expression.

* Correlation with the ground truth deseq results in order to find the similarities in gene expression.

* Analysing the results in suitable graphs and conclusions.


## 📃 Files in our system:

* Coverage.py - Coverage simulator that prepers the data, removes low quality data, chooses the random samples and perform the coverage generations.

* deseq.Rmd - Performs the deseq on the rawcounts and metadata and calculates the correlation with ground truth. Creates the table of all the different runs of sampling and coverage.

* graphs.py - Analzise the data in the table in graphs to allow us to understand the results and answer our questions.


## 🐭 The data:
The data was taken from the experiment "Aging increases cell-to-cell transcriptional variability upon immune stimulation" that was published in "Science" by the AAAS. Our data consists of RNA taken from mice's T cells and was organized by gene reads. Meaning, in order to analyze sequencing coverage, we'll use the genes reads instead of single nucleotides reads.


##	👥 Authors:
[Noa Kirsh](https://github.com/NoaKirsh) and Itai Friedman.
As a part of Bioinformatics project 236523 in the Technion.
