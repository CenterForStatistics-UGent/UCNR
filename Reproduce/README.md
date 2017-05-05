## Contents

This folder contains:
1. Data folder: The Neuroblastoma dataset underlying the simulation study.
2. Scripts folder: Some auxiliary scripts that allow to read in the Neuroblastoma dataset and functions to obtain plots as presented in the manuscript.
3. RobustnessStudy scripts: These scripts can be run to obtain simulation results as described in the paper, using either a Wilcoxon test or a t-test. 
4. Figures script: After the robustness simulation results have been obtained, those results can be used to generate plots such as those presented in the manuscript.

## How to obtain additional results

The manuscript and its appendices contain multiple simulation scenarios. The scripts within this folder can be altered easily to obtain these results, e.g.:
1. Scalability of the method: we performed simulations with 100 and 200 miRNAs (as opposed to 50 in the main manuscript simulation) to demonstrate the scalability of the method. These results may be obtained by changing the "testNrMis" variable accordingly.
2. Results for the different genes (as presented in supplementary material) may be obtained by altering the arguments of e.g. the plot.pvalues function in the Figures script.