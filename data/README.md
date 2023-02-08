# **Dataset**
The dataset can be found here <https://1drv.ms/u/s!AtMxC-h5vf1-iqYu3-xZMZaSjAkwFw?e=4LluyQ> (total size is ~12 Gb). The directory contains a folder for each tumor sample. Each tumor folder contains the following files:
- X_freebayes.vcf : Positions called as SNVs by freebayes algorithm.
- X_mutect.vcf : Positions called as SNVs by mutect algorithm.
- X_vardict.vcf : Positions called as SNVs by vardict algorithm.
- X_varscan.vcf : Positions called as SNVs by varscan algorithm.
- X_truth.bed : The truth set of spiked-in SNVs. The file format uses three columns: chromosome, genomic start, genomic stop.

Read here on how <https://www.synapse.org/#!Synapse:syn312572/wiki/62018> the synthetic tumor data was generated.

The directory also contains two real tumor dataset (real1 and real2), where the ground truth is known (through manual curation). You have data on all mutations in the real1 sample. However, real2 is split into two parts: in part1 (real2_part1) you have both VCF data and truth labels (bed file); in part2 (real2_part2), you only have VCF data. 

The final challenge is to predict the truth labels in real2_part2 as accurately as possible. Accuracy will be quantified using the F1-score <https://en.wikipedia.org/wiki/F-score> (a script/markdown to calculate this score is available here: <https://1drv.ms/u/s!AtMxC-h5vf1-iqdLPzI3tQxBIx_61Q?e=2flGaD> ).

## Additional data and resources
If you have time, you may add any additional sources of data that can help you improve your model, for example:

- Repeats, listing repeat regions in the human genome 
  - info <http://www.repeatmasker.org> ,
  - data <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz> ,
  - specs <http://hgdownload-test.cse.ucsc.edu/goldenPath/cioSav1/database/rmsk.sql>
- Mappability, quantifying sequence uniqueness in the human genome <http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/>

Note: All data is based on the version/build of the human genome referred to as **GRCh37/hg19**. So any additional data you use must be based on this genome build.
