# FoB-Project- This project is an introduction to the basic theory and practice of solving common problems
in bioinformatics. It is a masters group project at VU Amsterdam. 

Nonsynonymous (missense) mutations occur where a single nucleotide in a DNA codon is
substituted for another (a single nucleotide polymorphism or SNP), resulting in a change to
the amino acid that the codon codes for. These missense mutations can have no or little
effect on protein function and phenotype, but can sometimes result in significant changes
that can cause disease. The impact prediction tools PolyPhen-2 (Adzhubei et al., 2010)
and SIFT (Vaser et al., 2016) are designed to predict which mutations in DNA will cause
changes in the cell. They can be used to help interpret mutation data from patients who have
genetic diseases, but these methods must be validated to assess how accurate their
predictions are. Validation of these tools requires experimental data with accurate
annotations (i.e., a benchmark or gold standard dataset) against which we can compare the
performance of our tools. In this case, this means SNP data with annotations to indicate
whether they are benign or pathogenic to compare against the impact predictions of SIFT
and PolyPhen. We will use the database ClinVar for our gold standard dataset. Once we
have benchmarked the predictions, we can visualise the performance of these tools by
creating a ROC (Receiver Operating Characteristic) plot, which plots the True Positive Rate
(TPR) against the False Positive Rate (FPR). 
