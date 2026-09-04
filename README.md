# Table of contents
1. [Introduction](#introduction)
2. [What's the paper about](#whats-the-paper-about)
    1. [Methods](#methods)
    2. [Key findings](#key-findings)
    3. [Takeaway](#takeaway)
3. [Theoric prerequisites](#theoric-prerequisites)
    1. [Compositional data and spurious correlation](#compositional-data-and-spurious-correlation)
    2. [CLR transformation](#clr-transformation)
    3. [Sparsity and pseudocounts](#sparsity-and-pseudocounts)
    4. [Correlation matrixes](#correlation-matrixes)
    5. [NorTA (Normal to Anything)](#norta-normal-to-anything)
    6. [Pielou diversity index](#pielou-diversity-index)
4. [Research goal](#research-goal) 
5. [Structure of the repository](#structure-of-the-repository)
    1. [Folder explanations](#folder-explanations)
    2. [Repository tree structure](#repository-tree-structure)
    3. [Scripts in 02_new_scripts](#Scripts_in_02_new_scripts)
6. [Versions](#versions)
7. [Syntax](#syntax)


# Introduction
 
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

The paper investigates the biases affecting correlation estimates in metagenomic data, which arise from the compositional nature of the data, within-sample diversity, and high sparsity (many unobserved/zero taxa).    Using simulated data, the authors show that standard compositional-data transformations, especially the centered log-ratio (CLR)   allow simple Pearson correlation to reliably recover the true correlation structure, particularly in high-dimensional settings.  
Sparsity, however, remains an open issue, tending to underestimate negative correlations.
 
### Methods
- Gaussian data simulated using the R package `mvtnorm`, varying dimensionality (D) and within-sample diversity (*Pielou* index, for more information about it see chapter 3.6, [Pielou diversity index](#pielou-diversity-index)), to isolate compositional biases under L1 vs. CLR normalization.
- Realistic sparse data simulated with the "Normal to Anything" (NorTA) approach, using a zero-inflated negative binomial marginal distribution ([zero-inflated negative binomial Wikipedia page](https://en.wikipedia.org/wiki/Zero-inflated_model)); zeroes replaced with 65% of the detection limit before CLR.
- Validation on real data from the HMP2 gut microbiome dataset (single healthy subject: *69-001*, contributed with 51 samples in total, the highest count in the dataset), comparing CLR + Pearson against Rho and SparCC.  

### Key findings
- L1-normalized correlations are strongly biased by within-sample diversity and do not improve with dimensionality.
- CLR-normalized correlations are independent of diversity, and their bias decreases rapidly with dimensionality (negligible above ~100 taxa), typical of metagenomic datasets.
- On real HMP2 data, CLR + Pearson closely matches Rho and SparCC at high dimensionality (OTU level), while differences grow at low dimensionality (phylum level).
- Data sparsity remains a limitation: error increases with the fraction of zeroes, especially for negative correlations; CLR mitigates but does not eliminate this bias.  

### Takeaway
In typical high-dimensional metagenomic settings, simple Pearson correlation on CLR-transformed data is a robust and computationally cheap choice, with more complex compositional-correction methods becoming relevant mainly at low dimensionality.  
 Sparsity handling remains the main open challenge.  

## Theoric prerequisites
This section covers briefly the theoric concepts needed to understand the whole project.

### Compositional data and spurious correlation

Compositional data are vectors of proportions that sum to 1 (by a mathematical point of view they lie on a simplex).  
Because of this constraint, standard Euclidean distances and correlations may not properly capture the relationships between components, potentially leading to **spurious correlations**.  
More information can be found on the [compositional data Wikipedia page](https://en.wikipedia.org/wiki/Compositional_data).


### CLR transformation

The Centered Log-Ratio (CLR) transformation maps compositional data out of the constrained simplex space into an unconstrained, real-valued space, allowing tools like Pearson correlation ([Pearson correlation Wikipedia page](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)) to be applied reliably.

$$CLR(x_i) = \ln\left( \frac{x_i}{g(x)} \right) = \ln(x_i) - \frac{1}{D}\sum_{j=1}^{D}\ln(x_j)$$

where 

$$g(x) = \left( \prod_{i=1}^{D} x_i \right)^{1/D} $$

is the geometric mean.  

The second form (log of each value minus the row-wise mean of the logs) is the one implemented throughout this repository, e.g. in `clr_pearson.R`:

```r
log_x <- log(x)
log_x - rowMeans(log_x)
```

Note: CLR requires $`x_i > 0`$, so zero counts must first be replaced with pseudocounts (see `pseudocount.R`).

### Sparsity and pseudocounts

One of the main problems of metagenomic data is their sparsity: the dataset has many zeroes.  
A zero doesn't mean an absence of microorganisms, it means the detection method used to collect the data couldn't find them (their abundance is below the detection limit).  
For a better understanding of this problem, here is a useful page: [Viable But Nonculturable (VBNC) - Wikipedia](https://en.wikipedia.org/wiki/Viable_but_nonculturable).  
Pseudocounts are small positive values used to replace zeros before applying a log-based transformation like CLR  (since `log(0)` is undefined).  
For reference : [Wikipedia page for Pseudocounts](https://en.wikipedia.org/wiki/Additive_smoothing#Pseudocount).  
Honorable mention for a better understanding of what a pseudocount is:
[more intuitive explanation from Reddit](https://www.reddit.com/r/bioinformatics/comments/8acfss/can_someone_give_a_pseudocounts_for_dummies/).


### Correlation matrixes

A correlation matrix is a table of numbers where each entry gives the correlation strength between two variables. It has the following properties:

1. **Symmetric**: the correlation between X and Y is the same as between Y and X (correlation is commutative).
2. **Values range from -1 to +1**: 0 means no linear correlation, +1/-1 mean perfect positive/negative correlation.
3. **Positive semi-definite (PSD)**: all eigenvalues must be ≥ 0.
4. **Diagonal values equal to 1**: a variable is perfectly correlated with itself.

Not every matrix satisfying properties 1, 2, and 4 is automatically PSD, see `generate_matrix_factors.R` for how this repository builds matrixes that are PSD by construction.  
For a more detailed reference: [what_is_a_correlation_matrix](https://www.displayr.com/what-is-a-correlation-matrix/)

### NorTA (Normal to Anything)

NorTA is a simulation technique used to generate correlated data with **any target marginal distribution**, starting from correlated Gaussian data.  
It works in two steps:

1. **Correlated normal data**: draw samples from a multivariate normal distribution with the desired correlation matrix (`rmvnorm()`).
2. **Marginal transformation**: convert each normal value to its percentile via the normal CDF (`pnorm()`), then map that percentile through the quantile function of the target distribution (e.g. `qzinegbin()` for sparse, zero-inflated counts).

Since both steps are **monotonic transformations**, rank correlation is preserved, but not exactly the Pearson correlation, since that is scale-sensitive and the mapping is non-linear. This is why the correlation recovered from the final simulated counts is not identical to the one used to generate the underlying normal data.

Implemented in `NorTa_simulation.R` (fixed zero-inflation shared by all taxa) and `data_sim_ph_driven.R` (per-taxon zero-inflation driven by an environmental pH gradient).

### Pielou diversity index

Pielou index, usually denoted as P(x) or J', describes how equally distributed the abundances of different species (or taxa) are within a given sample. It is calculated as:

$$P(x)= \frac{H(x)}{\ln(D)} $$

where D is the dimensionality of the sample (the number of taxa considered) and H(x) is the Shannon entropy index:

$$H(x) = -\sum_{i=1}^{D} p_i \cdot \ln(p_i)$$

**What it measures, and what it doesn't.**  
Pielou is an *evenness* index, not a *richness* index: it does not tell you how many species are present (that's what D itself, or richness metrics, are for), but rather how uniformly the abundance is spread across the species that are present.  
Two samples with the same D can have very different Pielou values: if every taxon has roughly the same abundance, P(x) is close to 1; if abundance is dominated by one or a few taxa while the rest are rare, P(x) drops toward 0, even though richness hasn't changed.  
Implemented in `pielou_ind.R`.  
The script `pielou_ind.R` is used in `ph_pielou_confrontation.R` to generate the bar plots.

## Research goal

This project reproduces the paper's core pipeline (filtering, CLR, Pearson correlation) on the same real data (HMP2, subject 69-001), validating it against a value reported in the paper (Pielou index ≈ 0.68).   
It also fixes a methodological gap found along the way: the standard method to simulate a correlation matrix with controlled sparsity (random values corrected with `nearPD()`) destroys nearly all the imposed zeroes, so this repository builds correlation matrixes that are valid by construction instead (`generate_matrix_factors.R`).  
Then, it proposes a new approach to sparsity, the paper's main open problem: instead of a fixed zero-inflation parameter, sparsity here depends on an environmental driver (pH) and each taxon's ecological niche (`data_sim_ph_driven.R`).  
Finally, the project allows to generate a small family of bar plots that
allow to observe how the Pielou Index varies in function of the (fake) taxa's ph tolerances.

## Structure of the repository

The structure of this repository was originally designed to replicate the original repository from the paper mentioned above.  
As the project evolved, I decided to keep some folders and separate the newly written scripts from the ones in the original repository.  
The folders named *"00_data"* and *"01_from_paper"* contain the original scripts, while *"02_new_scripts"* and *"03_discarded_methods"* contain the newly written scripts for this project.  
Why keep the old scripts from the original repository?  
The original work was excellent, and understanding it was a good exercise in itself.  
Very few changes were made, the most important being the use of the `here` R package, which makes the scripts more reproducible and less dependent on hard-coded paths.  
The rest of the work consisted of restyling the scripts according to the Tidyverse guidelines (see [Syntax](#syntax)) and writing test scripts for them.  


### Folder explanations

***Note about the scripts: every script has a header stating what the script does, what are its inputs and outputs and eventually few notes which covers what differs from the previous voices.***  

Here follows a brief explanation of every folder:

- **00_data**: contains the raw `.csv` data and the `dataprocess.R` script, which generates the `.rds` version of the data.   
The choice of "00" is conceptual: everything starts from data (yeah, I know it may not sound 100% true, but I thought it was clever).

- **01_from_paper**: contains all the original scripts.  
Their order mirrors the one in the original repository.

- **02_new_scripts**: contains all the new scripts that do not come from the original paper.   
The "demo" subfolder contains demonstrative scripts for some of the other scripts inside `02_new_scripts` (the idea was to write them so they could explain some of the most complicated scripts in that folder).

- **03_discarded_methods**: a somewhat misleading folder name, it contains a couple of scripts that I abandoned because they turned out not to be useful for the project as a whole.

- **outputs**: present in the original repository; stores the outputs coming mainly from `01_from_paper/` (like in the original project) which are not necessarily plots.

- **Plots**: same idea as the outputs folder, except it stores mainly the plots.

- **requirements**: contains one script listing all the libraries used in this repository (in the original repo it was just a loose script; I decided to put it in its own folder).  

*Note about the libraries: every library used, in every script, has a link to its CRAN/GitHub page and a brief description of what it does, written in two comment lines above the `library(<library-name>)` command.*

- **test**: contains all tests, organized into subfolders that mirror the structure of the folders containing the scripts under test.

### Repository tree structure
Here follows the tree structure:


```
.
├── 00_data/
│   ├── dataprocess.R
│   ├── meta_HMP2.rds
│   ├── otu_HMP2.rds
│   ├── README_DATA.md
│   ├── taxonomy.rds
│   └── raw/
│       ├── meta_HMP2.csv
│       ├── otu_HMP2_16S.csv
│       └── taxonomy_HMP2_16S.csv
├── 01_from_paper/
│   ├── biases_parameters_examples/
│   │   └── biases_parameters_examples.R
│   ├── compositional_effects/
│   │   ├── 00_compute_error.R
│   │   └── 01_plots_compositional_effects.R
│   ├── correlation_sparsity/
│   │   └── correlation_sparsity_script.R
│   ├── method_comparison/
│   │   ├── CLR.R
│   │   ├── layout_signed.R
│   │   ├── methods_comparison.R
│   │   ├── methods_comparison_review.R
│   │   └── TRIU.R
│   ├── NorTa_example/
│   │   └── norta_example_script.R
│   ├── sparsity_effects/
│       ├── trials/
│           ├── 0_compute_error_htrlnorm_randzero.R
│           ├── 0_compute_error_zinb.R
│           ├── 0_compute_error_zinb_rare.R
│           ├── intermediate_Plots_htrlnorm.R
│           └── intermediate_plots_htrlnorm_randzero.R
│       ├── a_example_toy.R
│       ├── B0_compute_error_htrlnorm.R
│       ├── B0_compute_error_zinb_rare.R
│       ├── B1_plots_htrlnorm.R
│       ├── B1_plots_zinb_rare.R
│       ├── C_example_zero_on_CLR.R
│       └── Z_final_plot.R        
│   ├── zero_replacements/
│       └── zero_replacements_script.R
├── 02_new_scripts/
│   ├── demo/
│   │   ├── demo_clr_pearson.R
│   │   ├── demo_data_sim_ph_driven.R    
│   │   ├── demo_datasummary.R
│   │   ├── demo_filt_data.R
│   │   ├── demo_generate_matrix_factors.R
│   │   ├── demo_NorTa_simulation.R
│   │   ├── demo_pielou_ind.R
│   │   └── demo_pseudocount.R
│   ├── clr_pearson.R
│   ├── compositional_bias_D_P.R
│   ├── data_sim_ph_driven.R
│   ├── datasummary.R
│   ├── filt_data.R
│   ├── generate_matrix_factors.R
│   ├── NorTa_simulation.R
│   ├── ph_pielou_confrontation.R
│   ├── pielou_ind.R
│   └── pseudocount.R
├── 03_discarded_methods/
│   ├── demo_clr_pearson_old.R
│   └── generate_matrix_with_zeroes.R 
├── outputs/
│   └── methods_comparison_outputs/
├── Plots/
├── requirements/
│   └── requirements.R
└── test/
    ├── test_00_data/
    │   └── test-data_process.R
    ├── test_01_from_paper/
    │   ├── test_biases_parameters_examples/
    │   │   └── test-biases_parameters_examples.R
    │   ├── test_compositional_effects/
    │   │   ├── test-00_compute_error.R
    │   │   └── test-01_plots_compositional_effects.R
    │   ├── test_correlation_sparsity/
    │   │   └── test-correlation_sparsity_script.R
    │   ├── test_method_comparison/
    │   │   ├── test-CLR.R
    │   │   ├── test-layout_signed.R
    │   │   ├── test-methods_comparison.R
    │   │   ├── test-methods_comparison_reviews.R
    │   │   └── test-TRIU.R
    │   ├── test_norta_example/
    │   │   └── test-norta_example_script.R
    │   ├── test_sparsity_effects/
    │   │   ├── test_trials/
    │   │   │   ├── test-0_compute_error_htrlnorm_randzero.R
    │   │   │   ├── test-0_compute_error_zinb.R
    │   │   │   ├── test-0_compute_error_zinb_rare.R
    │   │   │   ├── test-intermediate_Plots_htrlnorm.R
    │   │   │   └── test-intermediate_Plots_htrlnorm_randzero.R
    │   │   ├── test-a_example_toy.R
    │   │   ├── test-B0_compute_error_htlrnorm.R
    │   │   ├── test-B0_compute_error_zinb_rare_ecdf.R
    │   │   ├── test-B1_plots_htrlnorm.R
    │   │   ├── test-B1_plots_zinb_ecdf.R
    │   │   ├── test-C_example_zero_on_CLR.R
    │   │   └── test-final_plot.R
    │   └── test_zero_replacements/
    │       └── test-zero_replacements_script.R
    ├── test_02_new_scripts/
    │   ├── test-demos/
    │   │   ├── test-demo_clr_pearson.R
    │   │   ├── test-demo_data_sim_ph_driven.R
    │   │   ├── test-demo_datasummary.R
    │   │   ├── test-demo_filt_data.R
    │   │   └── test-demo_pielou_ind.R
    │   ├── test-clr_pearson.R
    │   ├── test-compositional_bias_D_P.R
    │   ├── test-data_sim_ph_driven.R
    │   ├── test-datasummary.R
    │   ├── test-filt_data.R
    │   ├── test-generate_matrix_factors.R
    │   ├── test-NorTa_simulation.R
    │   ├── test-ph_pielou_confrontation.R
    │   └── test-pseudocount.R
    └── test_03_discarded_methods/
        └── test-generate_matrix_with_zeroes.R  
```

*Note: the folders "outputs" and "Plots" contain files that are not listed in the tree structure (which are the plots and outputs of all the scripts of this project).  
This choice was deliberate, the aim was to reduce the visual impact of the tree structure (which is already pretty big).  
So, why not simply put the unlisted files in the .gitignore? Good question.  
Let's say that this is to carry on the will to mirror the original folder.*

## Scripts in 02_new_scripts

Ok, this is where the fun goes live.  
Most of these scripts define a single function and are meant to be sourced, not run on their own.  
The `demo/` subfolder shows how to actually use them (some demos are missing).  
`compositional_bias_D_P.R` is the only standalone analysis script in this folder (for the moment).  
We can say there are three main types of scripts: base functions, orchestrators, and standalone analyses.  
Here are their explanations:

**Base functions**
- `datasummary.R`: has the function `datasum()`, a quick summary of a dataset (dimensions, mean/median/sd, skewness, zero count, and determinant/min eigenvalue for square matrices).

- `filt_data.R`: has the function `filt_data()`, filters an OTU table by prevalence and abundance.

- `pseudocount.R`: has the function `pseudocount()`, replaces zero counts with row-specific pseudocounts.

- `generate_matrix_factors.R`: has the function `generate_matrix_factors()`, builds a valid correlation matrix (PSD by construction) with an exact, controlled zero pattern.

- `pielou_ind.R`: has the function `pielou_ind()`, computes the Pielou diversity index, works on real count data (zeroes included).

**Orchestrators** (combine the base functions above)
- `clr_pearson.R`: contains the function `clr_on_data()`, runs the full filtering → pseudocount → CLR → Pearson correlation pipeline.

- `NorTa_simulation.R`: brings the `norta_simulation()`, simulates a sparse dataset with a known correlation structure (NorTA approach), with a single zero-inflation probability shared by all taxa.

- `data_sim_ph_driven.R`: the key `data_sim_ph_driven()`, same NorTA approach, but zero-inflation is derived per-taxon from an environmental pH gradient instead of being fixed.

**Standalone analysis**
- `compositional_bias_D_P.R`: quantifies L1 vs CLR bias as a function of dimensionality (D) and diversity (P), replicating the paper's Analysis 1 on a reduced grid. Produces `compositional_bias_results.rds` and a heatmap.
- `ph_pielou_confrontation.R`: quantifies the variation of Pielou index in function of tolerances of the (artificially generated) taxa. Produces a 
combined plot named `ph_pielou_confrontation.png`.


**Demos** (`02_new_scripts/demo/`)  
Each demo shows one of the functions above in action, with explanatory `cat()` output and, where relevant, printed results.

## Versions

Programming language used was [R](https://www.r-project.org/), updated at the version 4.6.1.  
The Integrated Development Environment (IDE) used to write the scripts was [RStudio](https://posit.co/products/open-source/rstudio), updated at its version 2026.08.2 + 200.  
To write and edit the README.md was used [Visual Studio code](https://code.visualstudio.com/), updated at its latest version (1.136.0).

## Syntax
About the syntax, the conventions stated in the Tidyverse guide have been followed, they can be found here: 

[Tidyverse_guidelines](https://style.tidyverse.org/syntax.html)

A very useful tool (which has been used intensively) for styling the code syntax is the library ["styler"](https://styler.r-lib.org/).  
This library allows to style the active script in Rstudio, changing indentation and spaces of the code (by default it respects the Tidyverse conventions).  
Note about the tests: the library used to create the tests is "testthat".  
The conventions for writing the tests can be found here: 

[testthat_conventions](https://r-pkgs.org/testing-basics.html).

