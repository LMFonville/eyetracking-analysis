---
title: Tobii TX300 Eyegaze Data Preprocessing
author: Leon Fonville
date: 'December 05, 2018'
output: 
   html_document:
     keep_md: true
     toc: true
     toc_float: 
         collapsed: true
     code_folding: hide
     self_contained: true
---

# Introduction
This document reports a series of preprocessing steps used to classify eyegaze
data as fixations or saccades. The chosen processing steps are motivated by
selected reading from the literature, some software packages, and available
information on processing steps from Tobii. If I failed to cite any relevant
sources please let me know!

#Acknowledgements 



This document makes use of the dplyr, ggplot, and cowplot packages.
Many of the preprocessing steps are inspired by the Tobii documentation by 
[Olsen, 2012](https://www.tobiipro.com/siteassets/tobii-pro/learn-and-support/analyze/how-do-we-classify-eye-movements/tobii-pro-i-vt-fixation-filter.pdf/?v=2012)
The eyegaze velocity thresholding is based on work by [Mould et al., 2012](https://www.sciencedirect.com/science/article/pii/S0042698911004214) 
as well as the GazePath package by [van RensWoude et al., 2017](https://link.springer.com/article/10.3758%2Fs13428-017-0909-3).


```r
#load libraries
library(dplyr)
library(ggplot2)
library(cowplot)
```


# Summary


```r
#Some housekeeping for the summary table
nTrials = unique(input$trialIdx)
runIdx = unique(input$runIdx)
if (unique(input$clockCode==1)==TRUE){stimType="Analog Clock"
} else if (unique(input$clockCode)==2){stimType="Digital Clock"
}else{stimType="Invalid Code"}
if (unique(input$trialCode==1) | unique(input$trialCode==2)){trialType="Perception Trial"
} else if (unique(input$trialCode==3)){trialType="Imagery Trial"
} else{trialType="Invalid Code"}
#prep for reporting trials as subheadings, note the extra line before the closing quotation
subheader <- "### Run %d Trial number %d

"
```


 |                   |                        |
 |:------------------|-----------------------:|
 | Participant ID:   | 999 |
 | Experiment Run:   | 1             | 
 | Number of Trials: | 4    | 
 | Stimulus:         | Analog Clock           | 
 | Trial:            | Imagery Trial          | 
 
# Remove invalid data
We'll average the eyegaze position across both eyes. If no valid data is available for one eye, the data from the other eye will be used. If both eyes contain invalid data it will be excluded.






```r
#Calculate an average of eyegaze position across eyes when valid data is available
averageEye.Fun <- function(df){
  # Find left and right eye columns
  # Column index contains XY eyegaze data and validity code
  leftIdx <- grep("left", colnames(df))
  rightIdx <- grep("right", colnames(df))
  
  # Remove eyedata if valid code is greater than 1
  df[df[leftIdx[3]]>1, leftIdx[1:2]] = NA
  df[df[rightIdx[3]]>1, rightIdx[1:2]] = NA
  
  # Remove any eyedata outside of logical range (0-1)
  eyeIdx <- grep("Eye", colnames(df)) #columns that contain eyedata
  df[,eyeIdx] <- apply(df[,eyeIdx], 2, function(x) ifelse(x > 1 | x < 0, NA, x))
  
  # Remove "valid" data along one axis if other axis has invalid data
  NAidx <- which(is.na(df[,rightIdx[1]])!=is.na(df[,rightIdx[2]]))
  df[NAidx,rightIdx[1:2]] = NA
  NAidx <- which(is.na(df[,leftIdx[1]])!=is.na(df[,leftIdx[2]]))
  df[NAidx,leftIdx[1:2]] = NA
  df$idx = seq(1:dim(df)[1])
  
  # Remove any instances where the calculated distance from the eyetracker is implausible
  distances = sort(unique(df$distance))
  peakIdx <- which(diff(distances,1)>mean(diff(distances,1)))
  # This assumes there's very little variance overall with maybe a few peaks
  thresh = distances[peakIdx]
  naIdx <- which(df$distance %in% thresh)
  df$distance[naIdx] = NA
  # Plot validity codes for each eye for a quick overview of data quality
  plotTheme <- theme(line = element_blank(), legend.position = "none", 
                     axis.text = element_blank(), axis.ticks = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  #LOOK INTO THEME(ASPECT.RATIO) TO FIX COLOUR DISCREPANCIES
  lValid <- ggplot(df, aes(x = idx, y = ppID, fill=factor(leftValid))) + geom_tile() +
    scale_fill_manual(values = c("0"="green", "1" = "yellow", "2" = "red","3" = "red","4" = "red")) +
    labs(y = "left", x = "") + plotTheme #+ coord_fixed(ratio = 1000)
  rValid <- ggplot(df, aes(x = idx, y = ppID, fill=factor(rightValid))) + geom_tile() +
    scale_fill_manual(values = c("0"="green", "1" = "yellow", "2" = "red","3" = "red","4" = "red")) +
    labs(y = "right", x = "") + plotTheme #+ coord_fixed(ratio = 1000) 
  validPlot <- plot_grid(lValid, rValid, ncol = 1, nrow = 2) 
  plot(validPlot)
  cat('\n\n')
  
  # Calculate mean eyegaze position using valid data
  df <- df %>% mutate(X = rowMeans(df[c(leftIdx[1],rightIdx[1])], na.rm = T)) %>% 
    mutate(Y = rowMeans(df[c(leftIdx[2],rightIdx[2])], na.rm = T)) %>% 
    select(-leftIdx, -rightIdx)
  # Flip Y-axis so that the top right (X,Y) = (1,1)
  df$Y = 1 - df$Y
  return(df)
}
input <- averageEye.Fun(input)
```

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


Data contains 119 invalid datapoints (8.29%).

# Median Filtering
A sliding window with a windowlength of 3 frames was applied across the eyegaze data and each datapoint was replaced by the median value of those 3 frames.
The smoothed data is shown below in black, with the original data points in red.


```r
medianFilter.Fun <- function(df){
  # First, identify gaps in recordings and append eyegaze vectors accordingly
  # This is meant to identify missing frames and assumes there's not a lot of
  # different ranges of intervals between recorded frames
  val = unique(df$interval)
  # Which values are 'a lot' bigger than the smallest value 
  # (if most intervals are 120Hz find ones at around 80Hz or lower)
  jumps = which(round(val/min(na.omit(val)), digits=0)>1)
  appendDF = df
  colApp = grep("timestamp", colnames(df))
  
  # Ugly loop solution to append rows to even out intervals between frames
  for (i in 1:length(jumps)){
    rowIdx = which(appendDF$interval==val[jumps[i]])
    for (j in 1:length(rowIdx)){
      rowIdx = which(appendDF$interval==val[jumps[i]])
      nRows = round(appendDF$interval[rowIdx[j]]/min(na.omit(val)), digits = 0)-1
      appendDF = appendDF %>% add_row(ppID = rep(NA, times = nRows),.before = rowIdx[j])
      appendIdx = which(is.na(appendDF$ppID))
      for (h in 1:(colApp-1)){
        appendDF[appendIdx,h] = appendDF[rowIdx[j]-1,h]
      }
    }
  }
  df = appendDF
  
  # Fill missing values and smooth data using median filtering with a fixed windowsize of 3
  # Median function will return NA when values are missing so it won't fill large gaps
    #rewrite later as a function that takes a windowlength argument
  df$Xf <- df %>% group_by(trialIdx) %>% select(trialIdx, X) %>% mutate(Xlead = lead(X, 1)) %>% 
    mutate(X = X) %>% mutate(Xlag = lag(X, 1)) %>% ungroup() %>% select(-trialIdx) %>% 
    apply(., 1, function(x) median (x, na.rm=T))
  df$Yf <- df %>% group_by(trialIdx) %>% select(trialIdx, Y) %>% mutate(Ylead = lead(Y, 1)) %>% 
    mutate(Y = Y) %>% mutate(Ylag = lag(Y, 1)) %>% ungroup() %>% select(-trialIdx) %>% 
    apply(., 1, function(x) median (x, na.rm=T))
  df$Df <- df %>% group_by(trialIdx) %>% select(trialIdx, distance) %>% mutate(Dlead = lead(distance, 1)) %>% 
    mutate(D = distance) %>% mutate(Dlag = lag(distance, 1)) %>% ungroup() %>% select(-trialIdx) %>% 
    apply(., 1, function(x) median (x, na.rm=T))
  
  # Let's add a more intuitive timeframe variable and replace gaps with median values
  df$frameTime <- df %>% group_by(trialIdx) %>% select(trialIdx, timestamp) %>% mutate(Flead = lead(timestamp, 1)) %>% 
    mutate(F = timestamp) %>% mutate(Flag = lag(timestamp, 1)) %>% ungroup() %>% select(-trialIdx) %>% 
    apply(., 1, function(x) median (x, na.rm=T))
  df <- df %>% group_by(trialIdx) %>% mutate(frameTime = frameTime - first(frameTime))
  return(df)
}
input <- medianFilter.Fun(input)
```



### Run 1 Trial number 6

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotMedian-1.png)<!-- -->

### Run 1 Trial number 13

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotMedian-2.png)<!-- -->

### Run 1 Trial number 20

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotMedian-3.png)<!-- -->

### Run 1 Trial number 31

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotMedian-4.png)<!-- -->


# Velocity Thresholding
We will calculate the velocity at each frame as the distance between the preceding and subsequent datapoint converted into an angle. 
The Euclidean distance between points is calculated in millimeters and is then converted to a visual angle using the recorded distance from the screen. 
The velocity is reported in degrees per second. 


```r
angleV.Fun <- function(df, scr_x, scr_y){
  # Adjust the gaze points to reflect the widescreen used for recording
  # Calculate the Euclidean distance and velocity using the adjusted coordinates and
  # identify local maxima in angular velocity on each trial
  df <- df %>% mutate(Xmm = Xf * scr_x) %>% 
    mutate(Ymm = Yf * scr_y) %>% 
    mutate(Dmm = Df * 10) %>% 
    group_by(trialIdx) %>%
    # Euclidean distance between preceding and subsequent points
    mutate(ED = sqrt((lead(Xmm, 1) - lag(Xmm, 1)) ^ 2 +
                         (lead(Ymm, 1) - lag(Ymm, 1)) ^ 2)) %>%
    # Convert ED to a visual angle
    mutate(AV = ifelse(!is.na(ED) & (lag(Dmm, 1) > 0), (atan2(ED, Dmm) * 180 / pi) /
                         (lead(frameTime, 1) - lag(frameTime, 1)), NA)) %>%
    # Calculate Euclidean velocity as a sanity check, will show a linear relation to AV
    mutate(EV = sqrt((lead(Xmm, 1) - lag(Xmm, 1)) ^ 2 +
                       (lead(Ymm, 1) - lag(Ymm, 1)) ^ 2)
           / (lead(frameTime, 1) - lag(frameTime, 1))) %>%
    # Identify local maxima as points higher than preceding and subsequent points
    mutate(lMax = ifelse((AV > lag(AV, 1) & AV > lead(AV, 1)), 1, 0))
  return(df)
}
input <- angleV.Fun(input, scr_x, scr_y)
```

## Find Optimal Velocity Threshold
The next step is to find the speed threshold that optimally separates the distribution of saccades from the distribution of fixational eye movements and noise. 
We'll be using local regression models to find the optimum speed threshold and using cross validation to optimise the regression model.


```r
#Requires Angular Velocity (AV) and sampling rate (in Hz) as input
optThresh.Fun <- function(lmV, Hz){
  # Make dataframe of local maxima
  df <- data.frame(idx = 1:Hz,
  # Speed thresholds range from minimum to maximum recorded velocities in equidistant number of steps
                   thresh = seq(from = min(lmV), to = max(lmV), length.out = Hz),
                   threshCount = NA,
  # Null represents the frequency of local maxima exceeding the speed threshold across a uniform null distribution
                   nullCount = seq(from = length(lmV), to = 0, length.out = Hz),
                   gap = NA,
                   fit = NA)
  # Get the actual frequency of local maxima exceeding the threshold
  df$threshCount <- sapply(df$thresh, function(x) {length(which(lmV > x))})
  # Compute gap statistic
  df$gap <- df$nullCount - df$threshCount
  # Smooth the gap statistic using a local regression model with optimised smoothing parameter
  # Use k-fold cross validation to find the best smoothing parameter
  k <- 10
  folds <- sample(x = 1:k, size = length(lmV), replace = TRUE)
  hRange <- seq(from = 0.1, to = 0.9, by = 0.01)
  cvMtrx <- matrix(rep(x = NA, times = k * length(hRange)),
                   nrow = length(hRange), ncol = k)
  for(i in 1:length(hRange)) {
    for(j in 1:k) {
      modelFit <- loess(formula = gap ~ log(thresh), data = df[folds != j, ], span = hRange[i])
      modelPredict <- predict(object = modelFit, newdata = df[folds == j, ])
      cvMtrx[i, j] <- mean(abs(df$gap[folds == j] - modelPredict), na.rm = TRUE)
    }
  }
  cvAME <- rowMeans(cvMtrx, na.rm = TRUE)
  hPick <- which.min(cvAME)
  gapFit <- loess(formula = gap ~ log(thresh), data = df, span = hRange[hPick])
  df$fit <- predict(gapFit)
  # Check for multiple local maxima in optimal fit
  while (length(which(diff(diff(df$fit)>=0)<0)) > 1){
    cvAME = cvAME[-hPick]
    hPick <- which.min(cvAME)
    gapFit <- loess(formula = gap ~ log(thresh), data = df, span = hRange[hPick])
    df$fit <- predict(gapFit)
  }
  
  # Plot optimal model if need be
  plotDF <- data.frame(x = c(df$thresh, rev(df$thresh)),
                       y = c(df$threshCount, rep(0, Hz)),
                       z = rep(df$nullCount, 2))
  optFit <- ggplot() +
    # Add null line
    geom_line(data = df, aes(x = thresh, y = nullCount), linetype = "dashed", size = 1) +
    # Fill frequency local maxima as shape
    geom_polygon(data = plotDF, aes(x = x, y = y), alpha = 0.4) +
    # Add model fit
    geom_line(data = df, aes(x = thresh, y = fit, color = "red"), size = 1.5) +
    # Add gap line
    geom_line(data = df, aes(x = thresh, y = gap), color = "blue", linetype = "dotted", size = 1) +
    # Add optimal threshold
    geom_vline(xintercept = df$thresh[which.max(df$fit)], size = 2) + 
    theme_minimal() +
    scale_x_continuous(limits = c(0, max(df$thresh)), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = "Speed Threshold", y = "Frequency of Local Maxima") +
    theme(legend.position = "none") + ggtitle("Distribution of Local Maxima in Eyegaze Speeds")
  
  # Return the data a list
  opt = list()
  opt$cvMatrix <- cvMtrx
  opt$df <- df
  opt$value <- df$thresh[which.max(df$fit)]
  opt$h <- hRange[hPick]
  opt$plot <- optFit
  return(opt)
}
optParam <- input %>% filter(lMax==1) %>% pull(AV) %>% optThresh.Fun(., 120)
optParam$plot #plot the data
```

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/fitModel-1.png)<!-- -->


The observed frequency of local maxima exceeding the speed threshold is depicted in grey. 
The dashed line depicts the null distribution and the dotted blue line represents the gap between the null and the observed distributions.
The red line is a smoothed fit of the gap between the two and the black bar represents the optimal speed threshold for
distinguishing fixations from saccades. 

The optimal speed threshold for distinguishing fixations from saccades was 22.47

## Optimal Velocity Threshold

### Run 1 Trial number 6

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotThresh-1.png)<!-- -->

### Run 1 Trial number 13

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotThresh-2.png)<!-- -->

### Run 1 Trial number 20

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotThresh-3.png)<!-- -->

### Run 1 Trial number 31

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotThresh-4.png)<!-- -->


# Eyegaze Classification

### Run 1 Trial number 6

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotCluster-1.png)<!-- -->

### Run 1 Trial number 13

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotCluster-2.png)<!-- -->

### Run 1 Trial number 20

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotCluster-3.png)<!-- -->

### Run 1 Trial number 31

![](C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/reports/pp_999/pp999-run1-stim3-preproc_files/figure-html/plotCluster-4.png)<!-- -->


---
title: "preproc.R"
author: "leonf"
date: "Wed Dec 05 14:06:57 2018"
---
