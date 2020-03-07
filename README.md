# Bioinfornatics-project - 236523

Studying the influence of sequencing coverage / replicates on differential analysis.

Main questions:

In our project we will attempt to answer the following questions:
*	Primary: What is the ideal balance between coverage and replicates for a given budget?

*	Ideal balance depending on amount of budget

* Replicates Origin: Tradeoff between number of individuals

and number of replicates from each one.

## Work Flow:

* Load the data.

* Removing all the low quality data based on the correlation to the ERCC genes.

* Choosing random samples for a given number of mice and percentage of samples.

* Coverage - Generating random gene samplings from the data, based on their probability to be chosen in the “real world” for different.

* Creating rawcounts and metadata files for deseq, including all the genes and samples for a  given number of mice and coverage percentage.

* DESEQ - running deseq to analyse sequencing data and test for differential Expression.

* Correlation with the ground truth deseq results in order to find the similarities in gene expression.

* Analysing the results in suitable graphs and conclusions.
