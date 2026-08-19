# Table of contents
1. [Introduction: correlation-measure-replication](#introduction)
2. [What's the paper about](#paragraph1)
    1. [Methods](#subparagraph1)
    2. [Key findings](#subparagraph2)
    3. [Takeaway](#subparagraph3)
3. [Structure of the repository](#paragraph2)
4. [Research goal](#paragraph3)

# Introduction: correlation-measure-replication
 
The aim of this repository is to store the code needed to replicate and
extend the analysis from the paper *"Correlation Measures in Metagenomic
Data: The Blessing of Dimensionality"* (from now on I will refer to it as "the paper"):
 
> Fuschi, A.; Merlotti, A.; Tran, T.D.B.; Nguyen, H.; Weinstock, G.M.;
> Remondini, D. Correlation Measures in Metagenomic Data: The Blessing of
> Dimensionality. *Appl. Sci.* **2025**, *15*, 8602.
> https://doi.org/10.3390/app15158602

This repository is used for the final projects of both the **Software and
Computing for Applied Physics** and the **Statistical Data Analysis for Applied Physics** courses.

## What's the paper about

The paper investigates the biases affecting correlation estimates in metagenomic data, which arise from the compositional nature of the data, within-sample diversity, and high sparsity (many unobserved/zero taxa). Using simulated data, the authors show that standard compositional-data transformations — especially the centered log-ratio (CLR) — allow simple Pearson correlation to reliably recover the true correlation structure, particularly in high-dimensional settings. Sparsity, however, remains an open issue, tending to underestimate negative correlations.
 
**Methods:**
- Gaussian data simulated using the R package `mvtnorm`, varying dimensionality (D) and within-sample diversity (*Pielou* index, P, for more informations about it: [Pielou_index_description_link](https://www.statology.org/how-to-calculate-interpret-pielous-evenness-index/)), to isolate compositional biases under L1 vs. CLR normalization.
- Realistic sparse data simulated with the "Normal to Anything" (NorTA) approach, using a zero-inflated negative binomial marginal distribution; zeros replaced with 65% of the detection limit before CLR.
- Validation on real data from the HMP2 gut microbiome dataset (single healthy subject: *69-001*, contributed with 51 samples in total, the highest count in the dataset), comparing CLR + Pearson against Rho and SparCC.  

**Key findings:**
- L1-normalized correlations are strongly biased by within-sample diversity and do not improve with dimensionality.
- CLR-normalized correlations are independent of diversity, and their bias decreases rapidly with dimensionality (negligible above ~100 taxa) — typical of metagenomic datasets.
- On real HMP2 data, CLR + Pearson closely matches Rho and SparCC at high dimensionality (OTU level), while differences grow at low dimensionality (phylum level).
- Data sparsity remains a limitation: error increases with the fraction of zeros, especially for negative correlations; CLR mitigates but does not eliminate this bias.  

**Takeaway:** in typical high-dimensional metagenomic settings, simple Pearson correlation on CLR-transformed data is a robust and computationally cheap choice, with more complex compositional-correction methods becoming relevant mainly at low dimensionality. Sparsity handling remains the main open challenge.  

## Structure of the repository



## Research goal

This project builds on the simulation framework and findings described in the paper.  
The aim of this repository is to reproduce the results obtained in the paper and to add new proposal of solutions for the (persistent) problem of sparsity.  


# Versions

Programming language used was R, updated at the version 4.6.1.
The Integrated Development Environment (IDE) used to write the scripts was RStudio, updated at its version 2026.08.1 + 195.



# Syntax
About the syntax, the conventions stated in the Tidyverse guide have been followed, they are found here: 

https://style.tidyverse.org/syntax.html

A very useful tool (which has been used intensively) for styling the code syntax is the library ["styler"](https://styler.r-lib.org/).