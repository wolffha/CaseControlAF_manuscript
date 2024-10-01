# CaseControlAF_manuscript

This repository contains the code used to perform simulations to test the R package CaseControlAF (<https://github.com/wolffha/CaseControlAF/>)

The files are as follows:

* **simulation_source.R**: contains the source code (as a function) to generate simulated datasets
* **AFAnalysis_server.R**: contains the code to analyze and plot the results using CaseControlAF and the simulated dataset
* **AFSimulations.R**: contains the code used to run all simulations included in manuscript
* **CaseControlSE_correction_final.Rmd**: R notebook (with code and output) from development of CaseControl_SE correction framework
* **CaseControlSE_correction_final.html**: knitted notebook output as HTML from development of CaseControl_SE correction framework
* **CaseControlAF_ccGWAS.Rmd**: R notebook (with code chunks and plain text) for example usage of package methods with All of Us data to perform case-case GWAS
* **CaseControlAF_sexChromosomeSim.R**: Code to perform simulations to generate AFs for autosomes and sex chromosomes 
* **RootSimulations.R**: Code to perform simulations and output all root solution options for CaseControl_AF and plot results - toward showing only one root (x2) is within [0,1]

Please contact Hayley Wolff (hayley.wolff@cuanschutz.edu) for additional details or data
