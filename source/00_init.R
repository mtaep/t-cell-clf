# ------------------------------------------------------------------------------
# 36625 Algorithms in Bioinformatics - DTU Bioinformatics
# June 2018 - Leon Eyrich Jessen - jessen@bioinformatics.dtu.dk
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# This script will install the required packages for exercises in Deep Learning
# for the 3-week course "36625 Algorithms in bioinformatics" June 2018 at
# Department of Bio and Health Informatics, Technical University of Denmark
# ------------------------------------------------------------------------------



# Install required packages
# ------------------------------------------------------------------------------
install.packages("tidyverse")
install.packages("devtools")
devtools::install_github("omarwagih/ggseqlogo")
devtools::install_github("rstudio/keras")
library("keras")
install_keras()
