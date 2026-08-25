# List of all the packages necessary to reproduce the analysis with the associated
# command to install them.
# Important: this script is not intended to be run, it is intended as a list of
# packages used in the whole repository


# -------- List of packages --------

# cowplot: draw ggplot2 in new figures
# https://cran.r-project.org/web/packages/cowplot/index.html
install.packages("cowplot")
# doParallel: parallel backend for foreach
# https://cran.r-project.org/package=doParallel
install.packages("doParallel")
# doSNOW: parallel backend for foreach, supports progress bars via snow clusters
# https://cran.r-project.org/web/packages/doSNOW/index.html
install.packages("doSNOW")
# extrafont: used to use fonts other than the standard PostScript fonts
# https://cran.r-project.org/web/packages/extrafont/index.html
install.packages("extrafont")
# foreach: parallel foreach loops
# https://cran.r-project.org/package=foreach
install.packages("foreach")
# ggplotify: converts base R plots and grobs into ggplot2 objects (as.grob)
# https://cran.r-project.org/web/packages/ggplotify/index.html
install.packages("ggplotify")
# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
install.packages("ggpubr")
# ggVennDiagram: Venn diagram based on ggplot2
# https://cran.r-project.org/web/packages/ggVennDiagram/index.html
install.packages("ggVennDiagram")
# grid: low-level graphics, viewport and grob management
# https://cran.r-project.org/web/packages/grid/index.html
install.packages("grid")
# gridExtra: arrange multiple grid-based plots on a page
# https://cran.r-project.org/web/packages/gridExtra/index.html
install.packages("gridExtra")
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
install.packages("here")
# igraph: used to create and analyze network graphs
# https://cran.r-project.org/web/packages/igraph/index.html
install.packages("igraph")
# MASS: statistical modeling and analysis
# https://cran.r-project.org/package=MASS
install.packages("MASS")
# magick: advanced image processing and manipulation (read, write, transform images)
# https://cran.r-project.org/web/packages/magick/index.html
install.packages("magick")
# Matrix: tools for working with dense and sparse matrices, including nearPD()
# https://cran.r-project.org/web/packages/Matrix/index.html
install.packages("Matrix")
# mvtnorm: generates multivariate normal and t distributions
# https://cran.r-project.org/web/packages/mvtnorm/index.html
install.packages("mvtnorm")
# progress: displays text progress bars for long-running loops
# https://cran.r-project.org/web/packages/progress/index.html
install.packages("progress")
# qualpalr: generate distinct qualitative color palette
# https://cran.r-project.org/web/packages/qualpalr/index.html
install.packages("qualpalr")
# styler: automatically formats R code according to a consistent style guide
# https://cran.r-project.org/web/packages/styler/index.html
install.packages("styler")
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
install.packages("testthat")
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
install.packages("tidyverse")
# vegan: used to elaborates shannon entropy
# https://cran.r-project.org/web/packages/vegan/index.html
install.packages("vegan")
# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
install.packages("VGAM")
# zCompositions: imputation of zeros, left-censored and missing values in compositional data
# https://cran.r-project.org/web/packages/zCompositions/index.html
install.packages("zCompositions")


# -------- List of packages from GitHub --------

# propr: package with proportionality rho method
# https://github.com/tpq/propr
devtools::install_github("https://github.com/tpq/propr")
# pulsar: parallel utilities for lambda selection along a regularization path
# https://github.com/zdk123/pulsar
devtools::install_github("zdk123/pulsar")
# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
devtools::install_github("zdk123/SpiecEasi")
# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
devtools::install_github("https://github.com/Fuschi/ToyModel")
