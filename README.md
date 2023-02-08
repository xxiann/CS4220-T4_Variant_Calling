# CS4220-Team-4

## Task
You are working on a project to call somatic single nucleotide variants (SNVs) in cancer. You are running 4 different SNV calling algorithms that sometimes give different results. The four SNV algorithms are: 
- MuTect <https://www.broadinstitute.org/cancer/cga/mutect>, 
- VarScan <http://dkoboldt.github.io/varscan/>, 
- FreeBayes <https://github.com/ekg/freebayes>,
- VarDict <https://github.com/AstraZeneca-NGS/VarDictJava>. 

Your task is to construct a robust meta-approach/method that integrates the calls from the different methods in a way that achieves highest accuracy on an independent test dataset.

## Final test data
This dataset (real2_part2) consists of a set of variant calls for the 4 algorithms on real tumor data. You must compute your predictions [0 = no SNV, 1 = SNV] for each of the candidate SNVs (union of predicted SNVs from all 4 methods). Please submit your predictions of true mutations in the test data as a bed file together with your report.

## Recommended reading
- Two somatic SNV meta-callers:
  - <http://genomebiology.com/2015/16/1/197>
  - <https://academic.oup.com/bioinformatics/article/35/17/3157/5288515> 
Read these papers to get inspirations of new features or approaches you could add to your model.
Paper describing the DREAM challenge on SNV prediction, including a description of the synthetic datasets used in this project: <http://www.nature.com/nmeth/journal/v12/n7/abs/nmeth.3407.html>
