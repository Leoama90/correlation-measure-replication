# correlation-measure-replication
 
The aim of this repository is to store the code needed to replicate and
extend the analysis from the paper *"Correlation Measures in Metagenomic
Data: The Blessing of Dimensionality"* (from now on I will refer to it as "the paper"):
 
> Fuschi, A.; Merlotti, A.; Tran, T.D.B.; Nguyen, H.; Weinstock, G.M.;
> Remondini, D. Correlation Measures in Metagenomic Data: The Blessing of
> Dimensionality. *Appl. Sci.* **2025**, *15*, 8602.
> https://doi.org/10.3390/app15158602
 
This repository is used for the final projects of both the **Software and
Computing for Applied Physics** and the **Statistical Data Analysis for Applied Physics** courses.
 
## Research goal

The aim of the paper is to investigate the __biases__ affecting correlation measures used 
to reconstruct microbial networks from metagenomic data.
This project investigates how data sparsity (i.e., the fraction of
zero-valued entries) affects the accuracy of correlation estimates obtained
through the centered log-ratio (CLR) transformation combined with Pearson's
correlation, on simulated compositional data mimicking metagenomic
abundance matrices. The work builds on the simulation framework and
findings described in the paper.

# Versions

The R version used was the  4.5.3.
Rstudio version was 2026.06.0+242.

# Syntax
About the syntax I follow the conventions stated in the Tidyverse guide, found here: 

https://style.tidyverse.org/syntax.html
