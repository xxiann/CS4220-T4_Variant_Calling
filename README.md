# CS4220-Team-4

## Content of repo
- `function.R` = Contains F1 and preprocessing functions used
- `model.R` = Contains code used to train the final model
- `final_model.rds` = Contains the final trained model
- `predict.R` = File to generate predictions for real2_part1
- `parse-vcf-snv-wrapper.R` = Parser function

## To reproduce the model and prediction 
1. Obtain all the data files for `real1`, `real2_part1` and `real2_part2`

2. Run `model.R`, which parses and preprocess training data, and trains and saves the final model 
- (change file names accordingly)

3. Then run `predict.R`, which loads the final model, parses and preprocess test data to obtain the final predictions.

## Task
You are working on a project to call somatic single nucleotide variants (SNVs) in cancer. You are running 4 different SNV calling algorithms that sometimes give different results. The four SNV algorithms are: 
- MuTect <https://www.broadinstitute.org/cancer/cga/mutect>, 
- VarScan <http://dkoboldt.github.io/varscan/>, 
- FreeBayes <https://github.com/ekg/freebayes>,
- VarDict <https://github.com/AstraZeneca-NGS/VarDictJava>. 

Your task is to construct a robust meta-approach/method that integrates the calls from the different methods in a way that achieves highest accuracy on an independent test dataset.

