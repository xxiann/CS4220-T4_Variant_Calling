# CS4220-Team-4

## Task
You are working on a project to call somatic single nucleotide variants (SNVs) in cancer. You are running 4 different SNV calling algorithms that sometimes give different results. The four SNV algorithms are: 
- MuTect <https://www.broadinstitute.org/cancer/cga/mutect>, 
- VarScan <http://dkoboldt.github.io/varscan/>, 
- FreeBayes <https://github.com/ekg/freebayes>,
- VarDict <https://github.com/AstraZeneca-NGS/VarDictJava>. 

Your task is to construct a robust meta-approach/method that integrates the calls from the different methods in a way that achieves highest accuracy on an independent test dataset.

## Content of repo
- function.R = Contains F1 and preprocessing functions used
- model.R = Contains code used for the model
- final_model.rds = Contains the trained model
- predict.R = File to generate predictions


