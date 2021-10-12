1_MuDiCA
================

## Step 1: Run Multiblock Discriminant Correspondence Analysis

``` r
# Set up block normalization
groups.norm.nom.bary.sum <- rowSums(groups.norm.nom.bary)
resDiCA.data <- (data.all.nom*matrix(data =groups.norm.nom.bary.sum, nrow = dim(data.all.nom)[1], ncol = length(groups.norm.nom.bary.sum),byrow=TRUE))

# Input data matrix for multiblock DiCA
rawData <- resDiCA.data

#_____________________________________________________________________
# Computations ----
# Run DiCA  ----

# Recode as factors ----
XYmat <- rawData
design.choose <- descriptors$AgeGen

g.masses <- NULL
# If you want to give equal weights to call the categories, use the g.masses code below
# g.masses <-  rep(1 / ncol(makeNominalData(XYmat)), length(unique(descriptors$AgeGen)))

resDiCA <- tepDICA(XYmat, make_data_nominal = FALSE, 
                   group.masses = g.masses,
                   #weight = rep(1, nrow(XYmat)),# -- if equal weights for all columns,                    
                   DESIGN = design.choose, graphs = FALSE)

# Inferences ----
set.seed(70301) # set the seed

# For random effects model so that we all have the same results.
nIter = 1000
resDiCA.inf <- tepDICA.inference.battery(XYmat,make_data_nominal = FALSE,
                                         DESIGN = design.choose,
                                         group.masses = g.masses,
                                         test.iters = nIter,
                                         #weight = rep(1, nrow(XYmat)), # -- if equal weights for all columns,
                                         graphs = FALSE)
```

    ## [1] "It is estimated that your iterations will take 0.34 minutes."
    ## [1] "R is not in interactive() mode. Resample-based tests will be conducted. Please take note of the progress bar."
    ## ================================================================================

### Save the output

``` r
save.image(file = "../data/fromMuDiCA.RData")
```
