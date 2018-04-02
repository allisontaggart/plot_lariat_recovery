# plot_lariat_recovery

Plots the steady-state lariat levels (defined as lariat reads mapped per
million hg19 reads) for all datasets.  Includes patient fibroblast RNAseq
data from all used genotypes (Controls, DBR1, TLR3, STAT1) with the all 
treatments (no treatment, ifna2b stimulation, pic stimulation, 24hrs HSV1
infection).  Plots each RNAseq experiment as a data point (each 
experiment is an average of 3 technical replicates) and calculates error
bars using SEM.

(Lariat reads mapped using my lariat read aligner:
https://github.com/allisontaggart/lariat, hg19 reads mapped using STAR w/
default parameters, and counting done with in-house scripts)

Full data/analysis as part of :
"Inborn Errors of RNA Lariat Metabolism in Humans with Brainstem Viral
Infection" - Cell. 2018 Feb 22;172(5):952-965.e18. doi: 10.1016/j.cell.2018.02.019.
