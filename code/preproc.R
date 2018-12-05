#' ---
#' title: Tobii TX300 Eyegaze Data Preprocessing
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
#' # Introduction
#'This document reports a series of preprocessing steps used to classify eyegaze
#'data as fixations or saccades. The chosen processing steps are motivated by
#'selected reading from the literature, some software packages, and available
#'information on processing steps from Tobii. If I failed to cite any relevant
#'sources please let me know!
#'
#' #Acknowledgements 
#+loadLib, echo=FALSE, message=FALSE, warning=FALSE
#' This document makes use of the dplyr, ggplot, and cowplot packages.
#' Many of the preprocessing steps are inspired by the Tobii documentation by 
#' [Olsen, 2012](https://www.tobiipro.com/siteassets/tobii-pro/learn-and-support/analyze/how-do-we-classify-eye-movements/tobii-pro-i-vt-fixation-filter.pdf/?v=2012)
#' The eyegaze velocity thresholding is based on work by [Mould et al., 2012](https://www.sciencedirect.com/science/article/pii/S0042698911004214) 
#' as well as the GazePath package by [van RensWoude et al., 2017](https://link.springer.com/article/10.3758%2Fs13428-017-0909-3).

#load libraries
library(dplyr)
library(ggplot2)
library(cowplot)
#'
#' # Summary
#+prepSum
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
#'
#'  |                   |                        |
#'  |:------------------|-----------------------:|
#'  | Participant ID:   | `r unique(input$ppID)` |
#'  | Experiment Run:   | `r runIdx`             | 
#'  | Number of Trials: | `r length(nTrials)`    | 
#'  | Stimulus:         | `r stimType`           | 
#'  | Trial:            | `r trialType`          | 
#'  

#' # Remove invalid data
#' We'll average the eyegaze position across both eyes. If no valid data is available for one eye, the data from the other eye will be used. If both eyes contain invalid data it will be excluded.
#+avgEye, results = "asis", dpi = 200, fig.height=1, fig.width=8,echo=FALSE, message=FALSE, warning=FALSE
#'
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
#'
#' Data contains `r (sum(is.na(input$X)) + sum(is.na(input$Y)))/2` invalid datapoints (`r round((((sum(is.na(input$X)) + sum(is.na(input$Y)))/2)/dim(input)[1])*100, digits = 2)`%).
#'
#' # Median Filtering
#' A sliding window with a windowlength of 3 frames was applied across the eyegaze data and each datapoint was replaced by the median value of those 3 frames.
#' The smoothed data is shown below in black, with the original data points in red.
#+smoothMedian, message=FALSE, warning=FALSE
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
#'
#+plotMedian, results = "asis", fig.width=8, fig.height=6,echo=FALSE, message=FALSE, warning=FALSE
for (i in 1:4) {
  # Overlay the filtered data on the raw data points for inspection
  current <- nTrials[i]
  cat(sprintf(subheader, runIdx, current))
  xPlot <- ggplot(input %>% filter(trialIdx==current)) + 
    geom_point(aes(x = frameTime, y = X, colour = "red"), show.legend=FALSE) +
    geom_line(aes(x = frameTime, y = Xf), size = 1, alpha=0.8) +
    labs(x = "", y = "X-axis") +
    theme_bw()
  yPlot <- ggplot(input %>% filter(trialIdx==current)) + 
    geom_point(aes(x = frameTime, y = Y, colour = "red"), show.legend=FALSE) +
    geom_line(aes(x = frameTime, y = Yf), size = 1, alpha=0.8) +
    labs(x = "Time elapsed", y = "Y-axis") +
    theme_bw()
  xyPlot <- plot_grid(xPlot, yPlot, ncol = 1, align = "v")
  plot(xyPlot)
  cat('\n\n')
}
#'
#' # Velocity Thresholding
#' We will calculate the velocity at each frame as the distance between the preceding and subsequent datapoint converted into an angle. 
#' The Euclidean distance between points is calculated in millimeters and is then converted to a visual angle using the recorded distance from the screen. 
#' The velocity is reported in degrees per second. 
#+getAngle
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
#' ## Find Optimal Velocity Threshold
#' The next step is to find the speed threshold that optimally separates the distribution of saccades from the distribution of fixational eye movements and noise. 
#' We'll be using local regression models to find the optimum speed threshold and using cross validation to optimise the regression model.
#+fitModel, fig.width=6, fig.height=4
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
#'
#' The observed frequency of local maxima exceeding the speed threshold is depicted in grey. 
#' The dashed line depicts the null distribution and the dotted blue line represents the gap between the null and the observed distributions.
#' The red line is a smoothed fit of the gap between the two and the black bar represents the optimal speed threshold for
#' distinguishing fixations from saccades. 
#' 
#' The optimal speed threshold for distinguishing fixations from saccades was `r round(optParam$value, digits = 2)`
#' 
#' ## Optimal Velocity Threshold
#+plotThresh, results = "asis", fig.width=8, fig.height=6,echo=FALSE, message=FALSE, warning=FALSE
# Label points as fixations or saccades using the optimal velocity threshold
input <- input %>% mutate(fixation = ifelse(AV < optParam$value, 1, 0)) 
# Temporarily remove NAs to classify clusters
naIdx <- which(is.na(input$AV))
input$fixation[naIdx] = 0
input <- input %>% mutate(cluster = ifelse(fixation==1, cumsum(1-fixation)+1,0))
input$fixation[naIdx] = NA
# Arrange the clusters into ascending values
input <- input %>% group_by(trialIdx) %>%
  mutate(cluster = plyr::mapvalues(cluster, unique(sort(cluster)), 0:(length(unique(cluster))-1)))
# Plot the thresholded timecourse
start = round(min(input$frameTime), digits = 0)
end = round(max(input$frameTime), digits = 0)
timeSeq = breaks = seq(from = start, to = end, by = end/10)

for (i in 1:4) {
  current <- nTrials[i]
  cat(sprintf(subheader, runIdx, current))
  xPlot <- ggplot(input %>% filter(trialIdx==current)) + 
    geom_line(aes(x = frameTime, y = Xf, 
                  colour = fixation, alpha = fixation), size = 1) +
    labs(x = "", y = "X-axis") +
    scale_x_continuous(breaks = timeSeq) +
    scale_alpha_continuous(range = c(1,1), na.value = 0) +
    theme_bw() + theme(legend.position = "none", 
                       panel.grid.minor.x = element_blank())
  yPlot <- ggplot(input %>% filter(trialIdx==current)) + 
    geom_line(aes(x = frameTime, y = Yf, colour = fixation, alpha = fixation), size = 1) +
    labs(x = "", y = "Y-axis") +
    scale_x_continuous(breaks = timeSeq) +
    scale_alpha_continuous(range = c(1,1), na.value = 0) +
    theme_bw() + theme(legend.position = "none")
  avPlot <- ggplot(input %>% filter(trialIdx==current)) + 
    geom_line(aes(x=frameTime, y = AV)) +
    labs(x = "Time elapsed", y = "Angular Velocity") + 
    scale_x_continuous(breaks = timeSeq) +
    geom_hline(yintercept = optParam$value, color = "red", linetype = "dashed") +
    theme_bw() + theme(legend.position = "none")
  xyvPlot <- plot_grid(xPlot, yPlot, avPlot, ncol = 1, align = "v")
  plot(xyvPlot)
cat('\n\n')
}
#'
#' # Eyegaze Classification
#+plotCluster, results = "asis", fig.width=8, fig.height=6,echo=FALSE, message=FALSE, warning=FALSE
# Filter the data to remove fixations shorter than 100ms
clusterL <- input %>% filter(cluster>0) %>% group_by(trialIdx, cluster) %>% summarise(rleL = n())
shortCL <- which(clusterL$rleL<(Hz/10))
shortIdx <- which(interaction(input$trialIdx, input$cluster) %in% interaction(clusterL$trialIdx[shortCL], clusterL$cluster[shortCL]))
input$cluster[shortIdx] = NA
clusterSeq = seq(from = start, to = end, by = end/(end*2))
# Plot the final fixations
for (i in 1:4) {
  current <- nTrials[i]
  cat(sprintf(subheader, runIdx, current))
  xyFSplot <- ggplot(input %>% filter(trialIdx==current & !is.na(fixation))) +
    geom_point(aes(x = Xf, y = Yf, size = AV, colour=frameTime),alpha = 0.5) +
    geom_path(aes(x = Xf, y = Yf),alpha = 0.1) +
    scale_x_continuous(breaks = seq(from=0, to = 0.9, by = 0.1), limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(from=0, to = 0.9, by = 0.1), limits = c(0,1), expand = c(0,0)) +
    facet_grid(fixation ~ .,
               labeller = labeller(fixation = c("0" = "saccades", "1" = "fixations"))) +
    labs(x = "X-Axis", y = "Y-axis") + theme_bw() + guides(size = FALSE) +
    scale_colour_gradientn(name = "time elapsed", breaks = clusterSeq, labels = as.character(clusterSeq) ,colours=rainbow(4)) +
    theme(legend.key.height = unit(0.5,'inch'))
  plot(xyFSplot)
  cat('\n\n')
}
