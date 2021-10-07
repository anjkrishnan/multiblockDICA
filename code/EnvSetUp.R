##############
## EnvSetUp ##
##############
## This is the code that sets up the environment (i.e., the very beginning 
## chunk in the original script). This script should be run at the beginning 
## of each script.
#_____________________________________________________________________
# Clean start ----
rm(list = ls())
graphics.off()
#_____________________________________________________________________
# parameters for pptx ----
leTitre   <- 'DICA: MHLA-c' # my title
leDir     <-  '../results/'       # Where am I
filename  <- 'GenAge-Block-DiCA'     # file name for results
path2save <-  paste0(leDir, filename)
#_____________________________________________________________________
#_____________________________________________________________________
# Preamble ----
leHome    <- getwd()
DataDir   <- '../data/'
# libraries ----

# install.packages('TInPosition') # if needed
# devtools::install_github('HerveAbdi/PTCA4CATA')
# devtools::install_github('HerveAbdi/data4PCCAR')

library(tidyverse)
library(ExPosition)
library(TExPosition)
library(TInPosition)
library(PTCA4CATA)
library(data4PCCAR)
library(Ckmeans.1d.dp)
library(ggrepel)

# A couple of functions we may need ----
function.path <- "functions/" # check
files.sources = list.files(function.path)
sapply(paste0(function.path, files.sources), source)
