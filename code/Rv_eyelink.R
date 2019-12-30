#' ---
#' title: Similarity of Fixation
#' author: Leon Fonville
#' date: '`r format(Sys.Date(), "%B %d, %Y")`'
#' output: 
#'    html_document:
#'      toc: true
#'      toc_float: 
#'          collapsed: true
#'      code_folding: hide
#'      self_contained: true
#' ---
#' 

#' In this document we will examine the overall similarity of eyegaze fixations
#' during perception and imagery trials. For each individual we will calculate
#' the overall density of fixations that fell within the ROI and calculate the
#' Rv coefficient for the density matrices. The Rv coefficient takes values
#' between 0 and 1 and is similar to the standard correlation coefficient.
#'
#' # Acknowledgements
#' Calculations are based on code in the FactoMineR package by 
#' [Le et al., 2008](https://www.jstatsoft.org/article/view/v025i01) 
#' and the equations in 
#' [Abdi, 2007](https://www.utdallas.edu/~herve/Abdi-RV2007-pretty.pdf).
#' This document makes further use of the dplyr, ggplot, reshape, and cowplot
#' packages.
#+loadLib, echo=FALSE, message=FALSE, warning=FALSE, include=FALSE

#Load libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

#' # Summary
#+prepSum
#Some housekeeping
ppID = unique(c(input1$ppID,input2$ppID))
if (length(unique(c(input1$runIdx, input2$runIdx))) > 1){
  runs = range(c(input1$runIdx, input2$runIdx))
  runIdx = paste0(runs[1], '-', runs[2])
} else {
  runIdx = unique(c(input1$runIdx, input2$runIdx))
}

nTrials <- paste0(length(unique(input1$trialIdx)), '-', length(unique(input2$trialIdx)))

if (unique(c(input1$clockCode, input2$clockCode)==1)==TRUE){stimType="Analog Clock"
} else if (unique(c(input1$clockCode, input2$clockCode)==2)==TRUE){stimType="Digital Clock"
} else{stimType="Invalid Code"}
if (exists('combo')){
  trialPair <- case_when(combo=='01' ~ 'FX_Stim',
            combo=='02' ~ 'FX_Img',
            combo=='12' ~ 'Stim_Img',
            combo=='03' ~ 'FX_Stim',
            combo=='04' ~ 'FX_Img',
            combo=='34' ~ 'Stim_Img')
}

Xcombo = strsplit(combo, "")[[1]][1]
Ycombo = strsplit(combo, "")[[1]][2]

if (Xcombo=='1' | Xcombo=='3'){
  Xtype = 'Stimulus'
} else if (Xcombo=='2' | Xcombo=='4'){
  Xtype = 'Imagery'
} else{Xtype = 'Fixation Cross'}

if (Ycombo=='1' | Ycombo=='3'){
  Ytype = 'Stimulus'
} else if (Ycombo=='2' | Ycombo=='4'){
  Ytype = 'Imagery'
} else{Ytype = 'Fixation Cross'}

#prep for reporting trials
subheader <- "### Run %d Trial number %d

"
#'  |                   |              |
#'  |:------------------|-------------:|
#'  | Participant ID:   | `r ppID`     |
#'  | Experiment Run:   | `r runIdx`   | 
#'  | Stimulus Type:    | `r stimType` | 
#'  | X Type Trials:    | `r Xtype`    |
#'  | Y Type Trials:    | `r Ytype`    |
#'  
# Get the limits of the ROI
if (!exists('imgGrid')){
  # size of the stimuli was 768*768 pixels and was fixed when presented
  imgGrid <- data.frame(x = c((res_x/2)-(768/2), (res_x/2)+(768/2), (res_x/2)+(768/2), (res_x/2)-(768/2), (res_x/2)-(768/2)),
                        y = c((res_y/2)-(768/2), (res_y/2)-(768/2), (res_y/2)+(768/2),(res_y/2)+(768/2), (res_y/2)-(768/2)))
}

#Get density for XY coordinates for each trial type
lim = c(unique(imgGrid$x), unique(imgGrid$y))

#Specify number of bins
n = 100

# Plot each matrix as a heatmap using ggplot and reshape
# Melt Var1 is the Y-axis as it collapses the matrix by columns
plotRaster.Fun <- function(m, t){
  ggplot(melt(m), aes(Var2 , Var1, fill= value)) +
  xlab("") + ylab("") + ggtitle(t) +
  geom_raster() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position = "none") + 
  scale_fill_viridis_c()
}

# Get the density for perception trials and imagery trials
X <- MASS::kde2d(x = na.omit(input1$Xf), y = na.omit(input1$Yf), n = n, lims = lim)$z # Perception
Y <- MASS::kde2d(x = na.omit(input2$Xf), y = na.omit(input2$Yf), n = n, lims = lim)$z # Imagery

# Store plots of raw density
densX.p <- plotRaster.Fun(X, "X")
densY.p <- plotRaster.Fun(Y, "Y")

# Scale the variables
Ys <- scale(Y, scale=FALSE)
Xs <- scale(X, scale=FALSE)

# Store plots of scaled density
dscaleX.p <- plotRaster.Fun(Xs, "X Scaled")
dscaleY.p <- plotRaster.Fun(Ys, "Y Scaled")

#' # Plot Density
#' Density was calculated within the ROI in a `r paste0(n,"-by-",n)` grid for
#' each trial type. `r Xtype` trials will be referred to as "X" and `r Ytype`
#' trials will be referred to as "Y"

#+plotDens, results='asis', fig.width = 12, fig.height = 8
plot_grid(densX.p, dscaleX.p, densY.p, dscaleY.p, ncol = 2, nrow = 2)

#' # Calculate Rv Coefficient
#' The Rv coefficient was defined as 
#' $$ Rv = \frac{trace(XX{^T}YY{^T})}{\sqrt{trace(XX^TXX^T)*trace(YY{^T}YY{^T})}} $$
#+plotMat, results='asis', fig.width = 12, fig.height = 4

# Create positive semi-definite matrices
XX <- Xs %*% t(Xs)
YY <- Ys %*% t(Ys)

#Store plots of matrices
mX.p <- plotRaster.Fun(XX, expression(XX^T))
mY.p <- plotRaster.Fun(YY, expression(YY^T))

plot_grid(mX.p, mY.p, ncol = 2)

# Plot the variance and covariance 
xyCov <- diag(x = diag(XX %*% YY), ncol = n, nrow = n)
xVar <- diag(x = diag(XX %*% XX), ncol = n, nrow = n)
yVar <- diag(x = diag(YY %*% YY), ncol = n, nrow = n)

xyCov.p <- plotRaster.Fun(xyCov, expression(XX^T~YY^T))
xVar.p <- plotRaster.Fun(xVar, expression(XX^T~XX^T))
yVar.p <- plotRaster.Fun(yVar, expression(YY^T~YY^T))

plot_grid(xyCov.p, xVar.p, yVar.p, ncol = 3)

#' # Output
#+Output
# Split the equation for easier notation
nom = sum(diag(XX %*% YY)) #scalar-valued covariance
denom = sqrt(sum(diag(XX %*% XX)) * sum(diag(YY %*% YY))) #sqrt of scalar-valued variances for X and Y
Rv = nom / denom

# Calculate the z criterion for a sampling distribution of the permutation coefficients

# Get the mean
xBeta <- (sum(diag(XX)))^2 / sum(diag(XX %*% XX))
yBeta <- (sum(diag(YY)))^2 / sum(diag(YY %*% YY))

# Get the mean of the set of permutated coefficients
E <- sqrt(xBeta * yBeta) / (n - 1) 

# Calculate the variance
xAlpha  <-  n - 1 - xBeta
yAlpha <-  n - 1 - yBeta
xDelta <- sum(diag(XX)^2) / sum(diag(XX %*% XX))
yDelta <- sum(diag(YY)^2) / sum(diag(YY %*% YY))
xC <- ((n - 1) * (n * (n + 1) * xDelta - (n - 1) * (xBeta + 2))) / (xAlpha*(n - 3)) 
yC <- ((n - 1) * (n * (n + 1) * yDelta - (n - 1) * (yBeta + 2))) / (yAlpha*(n - 3)) 
V <- xAlpha * yAlpha * ((2 * n * (n - 1) + (n - 3) * xC * yC) / (n * (n + 1) * (n - 2) * (n - 1)^3))

# Calculate z-score and accompanying one-sided(?) p-value
zscore = (Rv - E) / sqrt(V) 
pval = round(pnorm(-abs(zscore)), digits = 4)
if (pval==0){
  pval = "<0.001"
}

#collect the output as 1 row for a dataframe

if (exists('trialPair')){
  rvOut <- cbind(ppID, runIdx, stimType, trialPair, Rv, zscore, pval, nTrials)
} else {
  rvOut <- cbind(ppID, runIdx, stimType, Rv, zscore, pval)
}


#'  |                    |              |
#'  |:-------------------|-------------:|
#'  | Participant ID:    | `r ppID`     |
#'  | Experiment Run:    | `r runIdx`   | 
#'  | Stimulus Type:     | `r stimType` | 
#'  | Rv Coefficient:    | `r Rv`       | 
#'  | z-score:           | `r zscore`   | 
#'  | p-value:           | `r pval`     | 