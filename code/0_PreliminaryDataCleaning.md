0_Preliminary data cleaning
================

## Step 0 - Preliminary data cleaning

``` r
# Main dataset
data.master<- read.csv(paste0(DataDir, "MHLAc_AgeGen_AllData.csv"), header=TRUE)

# Rename some variables
# data.master$AgeGenCat <- gsub('25+.', '≥25', data.master$AgeGenCat)
# data.master$AgeGenCat <- gsub('18-24', '≤24', data.master$AgeGenCat)
# data.master$Age <- gsub('25+.', '≥25', data.master$Age)
# data.master$Age <- gsub('18-24', '≤24', data.master$Age)

# Set row names
rownames(data.master) <- data.master[,1]

# Specify categories
design.gender <- as.matrix(data.master[,4])
design.age <- as.matrix(data.master[,5])
design.age.gen <- as.matrix(data.master[,6])
design.age.sup <- as.matrix(data.master[,7])
design.age.gen.sup <- as.matrix(data.master[,8])
design.maj.clin <- as.matrix(data.master[,9])
design.clin <- as.matrix(data.master[,10])
design.major <- as.matrix(data.master[,11])

# Final data table
data.all <- as.matrix(data.master[,-c(1:11)])

# Nominalise data and design
data.all.nom <- as.matrix(makeNominalData(data.all))
design.gender.nom <- as.matrix(makeNominalData(design.gender))
design.age.nom <- as.matrix(makeNominalData(design.age))
design.age.sup.nom <- as.matrix(makeNominalData(design.age.sup))
design.age.gen.nom <- as.matrix(makeNominalData(design.age.gen))
design.age.gen.sup.nom <- as.matrix(makeNominalData(design.age.gen.sup))
design.major.nom <- as.matrix(makeNominalData(design.major))
design.clin.nom <- as.matrix(makeNominalData(design.clin))
design.maj.clin.nom <- as.matrix(makeNominalData(design.maj.clin))

# Get groups/blocks and make nominal
groups.norm<- read.csv(paste0(DataDir, "MHLAc_AgeGen_Blocks.csv"),header=TRUE)
groups.nominal <- makeNominalData(as.matrix(groups.norm[,1]))
colnames(groups.nominal) <-  sub(".", "", colnames(groups.nominal))

# Normalise data matrix by groups/tables
groups.norm.nom <- makeNominalData(as.matrix(groups.norm[,2]))
groups.norm.nom.bary <- scale(groups.norm.nom, center = FALSE, scale = colSums(groups.norm.nom))

# Organize the data and descriptors
descriptors <-  cbind(design.gender, design.age, design.age.gen, design.major, design.clin, design.maj.clin) 
colnames(descriptors) <- c("Gender", "Age", "AgeGen", "Major", "ClinCourse", "MajorClin")
descriptors <- data.frame(descriptors)
```

### Save the output

``` r
save.image(file = "../data/fromPrelim.RData")
```
