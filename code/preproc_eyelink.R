#' ---
#' title: Eyelink Data Preprocessing
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
#exampleSet = c(72, 91, 134, 136) # change later to a random sample
if ((input$trialCode[1]==1 | input$trialCode[1]==2)==TRUE){stimType="Analog Clock"
} else if ((input$trialCode[1]==3 | input$trialCode[1]==4)==TRUE){stimType="Digital Clock"
} else if ((input$trialCode[1]==0)==TRUE){stimType='Fixation Cross'
}else{stimType="Invalid Code"}
if (unique(input$trialCode==1) | unique(input$trialCode==3)){trialType="Perception Trial"
} else if (unique(input$trialCode==2 | unique(input$trialCode==4))){trialType="Imagery Trial"
} else if (unique(input$trialCode[1]==0)){trialType='Not Applicable'
} else{trialType="Invalid Code"}
#prep for reporting trials as subheadings, note the extra line before the closing quotation
subheader <- "### Run %d Trial number %d

"
#'
#'  |                   |                        |
#'  |:------------------|-----------------------:|
#'  | Tracker:          | "Eyelink"              |
#'  | Participant ID:   | `r unique(input$ppID)` |
#'  | Experiment Run:   | `r runIdx`             | 
#'  | Number of Trials: | `r length(nTrials)`    | 
#'  | Stimulus:         | `r stimType`           | 
#'  | Trial:            | `r trialType`          | 
#'  

#' # Remove invalid data
#' Any implausible eyegaze data will be removed and the Y-axis data will be inverted. 
#+validEye, message=FALSE, warning=FALSE
#'
# Calculate an average of eyegaze position across eyes when valid data is available

#THINK ABOUT IF ONLY ONE EYE IS AVAILABLE
validEye.Fun <- function(df, res_x, res_y){

  # Column index contains XY eyegaze data
  xIdx <- grep("X", colnames(df))
  yIdx <- grep("Y", colnames(df))
  tIdx <- grep("timestamp", colnames(df))
  # Remove any eyedata outside of logical range (1920 by 1280)
  df[,xIdx] <- sapply(df[,xIdx], function(x) ifelse(x > res_x | x < 0, NA, x))
  df[,yIdx] <- sapply(df[,yIdx], function(x) ifelse(x > res_y | x < 0, NA, x))
  
  # Remove "valid" data along one axis if other axis has invalid data
  NAidx <- which(is.na(df[,xIdx])!=is.na(df[,yIdx]))
  df[NAidx, c(xIdx,yIdx,tIdx)]  = NA
  
  # Flip Y-axis so that the top right (X,Y) = (1,1)
  df$X = df$eyeX
  df$Y = res_y - df$eyeY
  df <- df %>% select(-eyeX, -eyeY)
  
  # Downsample to a more sensible Hz to clean up the data
  # Trials were 3 seconds and will reduce Hz to 500Hz
  df2 <- df %>% group_by(trialIdx) %>% mutate(bin = ntile(timestamp, (Hz*3)/2)) %>%
    group_by(trialIdx, bin) %>% mutate_at(vars(X, Y, timestamp), mean) %>%
    distinct(X, Y, timestamp, .keep_all = T) %>% ungroup()
  
  return(df2)
}
input <- validEye.Fun(input, res_x, res_y)
#'
#' Data contains `r (sum(is.na(input$X)) + sum(is.na(input$Y)))/2` invalid datapoints (`r round((((sum(is.na(input$X)) + sum(is.na(input$Y)))/2)/dim(input)[1])*100, digits = 3)`%).
#' After removal of invalid datapoints all data was downsampled to `r Hz/2` Hz.
#'
#' # Median Filtering
#' A sliding window with a windowlength of `r win` was applied across the eyegaze data and each datapoint was replaced by the median value of those `r win` frames.
#' The smoothed data is shown below in black, with the original data points in red.
#+smoothMedian, message=FALSE, warning=FALSE
medianFilter.Fun <- function(df, win){
  
  # Find any jumps in data collection
  df <- df %>% group_by(trialIdx) %>% mutate(interval = timestamp - lag(timestamp,1)) %>% ungroup()
  #jumps <- which(scale(df$interval)>3)
  
  #val = unique(df$interval)
  # Which values are 'a lot' bigger than the smallest value 
  #jumps = which(df$interval > quantile(na.omit(val), 0.99))
  
  #df$jumps = 0
  #df$jumps[jumps] = 1
  
  # Identify gaps that should be ignored and gaps that can be filled
  # Tobii uses a max gap length of 75ms so we will stick with that
  maxGapL <- (Hz/2) * 0.075
  # Let's loop through trials to find gaps, it will take more time but we won't accidentally group NAs across trials.
  trials = length(nTrials)
  
  
  # We'll make a binary vector to indicate long gaps and add it the main dataframe
  gapfill = NULL
  for (t in 1:trials){
    subdf <- df %>% filter(trialIdx==nTrials[t])
    gaps = rle(is.na(subdf$X)) # Since we matched X and Y on NAs it doesn't matter which one we use
    gapdf = data.frame(idx = 1:length(gaps$values), 
                       trial = nTrials[t],
                       length = gaps$lengths,
                       values = gaps$values) 
    #add column indicating which gaps are too long to fill
    gapdf <- gapdf %>% mutate(omitGap = ifelse(values==TRUE & length >= maxGapL, NA, 1))
    # Append vector 
    for (r in 1:nrow(gapdf)){
      gapfill = c(gapfill, rep(gapdf$omitGap[r], gapdf$length[r]))
    }
  }
  df$gapfill = gapfill
  
  # Fill missing values and smooth data using median filtering with a variable windowsize
  # Median function will return NA when all values are missing so it won't fill large gaps
  
  windowSlide.Fun <- function(X, window){
    step = round((window-1)/2)
    
    # If we can't get a symmetrical window we will make the lag longer than the lead
    # So if the windowsize is 10 the value at point i will be the median of 
    if((window %% 2) == 0){
      Xlist <- c(lapply(1:step, function(x) lag(X, x)), lapply(1:(step-1), function(x) lead(X, x)))
    } else {Xlist <- c(lapply(1:step, function(x) lag(X, x)), lapply(1:step, function(x) lead(X, x)))}
    
    # Turn list into a dataframe
    Xdf <- do.call(data.frame, Xlist)
    
    # Add original timeseries
    Xdf$X <- X
    
    # Use median filter and return output
    Xf <- apply(Xdf, 1, function(x) median(x, na.rm = T))
    return(Xf)
  }
  
  # Apply filter to XY coordinates and timestamp
  df <- df %>% group_by(trialIdx) %>% 
    mutate(Xf = windowSlide.Fun(X, win), 
           Yf = windowSlide.Fun(Y, win), 
           frameTime = windowSlide.Fun(timestamp, win))
  df <- df %>% group_by(trialIdx) %>% mutate(frameTime = frameTime - first(frameTime),
                                             frameStep = frameTime-lag(frameTime,1)) %>% ungroup()
  # End by omitting the new values for the large gaps which should remain
  df <- df %>% mutate(Xf = ifelse(gapfill==1, Xf, Xf*gapfill),
                      Yf = ifelse(gapfill==1, Yf, Yf*gapfill ))
  return(df)
}
input <- medianFilter.Fun(input, win)
#'
#+plotMedian, results = "asis", fig.width=8, fig.height=6,echo=FALSE, message=FALSE, warning=FALSE
for (i in 1:length(nTrials)) { #length(nTrials)
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
#' Distance between pixels will be calculated and the Euclidean distance between points is then calculated
#'  in millimeters and is then converted to a visual angle using the recorded distance from the screen. 
#' The velocity is reported in degrees per second. Any extreme values (> 99.5th percentile) will be clipped as these most likely represent signal loss.
#+getAngle
angleV.Fun <- function(df, scr_x, scr_y, res_x, res_y, distance){
  
  # Convert pixel distance to angle
  df <- df %>% mutate(Xa = atan((scr_x / 2) / distance) * (180/pi)*2 / res_x * Xf,
                      Ya = atan((scr_y / 2) / distance) * (180/pi)*2 / res_y * Yf) %>%
  # Calculate velocity using distance between angles unless a jump occurred in this window
    # group_by(trialIdx) %>% mutate(AV = ifelse(jumps==0 & lead(jumps,1)==0 & lag(jumps,1)==0, 
    #                                           sqrt((lead(Xa,1) - lag(Xa,1))^2 + 
    #                                           (lead(Ya,1) - lag(Ya, 1))^2), NA)) %>%
    group_by(trialIdx) %>% mutate(AV = sqrt((lead(Xa,1) - lag(Xa,1))^2 +
                                              (lead(Ya,1) - lag(Ya, 1))^2)) %>%
  # Identify local maxima as points higher than preceding and subsequent points
    mutate(lMax = ifelse((AV > lag(AV, 1) & AV > lead(AV, 1)), 1, 0)) %>%
    ungroup() %>%
  # Clip any extremely large angle values as these most likely represent the tracker losing the eye
    mutate(AVf = ifelse(AV > quantile(na.omit(AV), 0.995), quantile(na.omit(AV), 0.995), AV))
  return(df)
}
input <- angleV.Fun(input, scr_x, scr_y, res_x, res_y, distance)
#' ## Find Optimal Velocity Threshold
#' The next step is to find the speed threshold that optimally separates the distribution of saccades from the distribution of fixational eye movements and noise. 
#' We'll be using local regression models to find the optimum speed threshold and using cross validation to optimise the regression model. 
#' Extreme outliers ( > 99.5th percentile) will be clipped before running optimisation.
#+fitModel, fig.width=6, fig.height=4, warning=FALSE
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
  cvAME_redux <- cvAME
  oldPick <- hPick
  # Check for multiple local maxima in optimal fit
  while (length(which(diff(diff(df$fit)>=0)<0)) > 1){
    # if exists, take the next smallest cvAME and accompanying h value
    cvAME_redux <-  cvAME[-oldPick]
    hPick <- which(cvAME %in% min(cvAME_redux))
    oldPick <- c(oldPick, hPick)
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
    scale_y_continuous(limits = c(0, max(df$nullCount)), expand = c(0,0)) + 
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
optParam <- input %>% filter(lMax==1) %>% pull(AVf) %>% optThresh.Fun(., Hz) #using clipped angular velocity
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
#' Any datapoints exceeding the speed threshold (red dotted line) are labelled as saccades and are depicted in dark blue. Fixations are shown in light blue, 
#'  but any fixations shorter than 100ms are discarded and are labelled with an intermediate shade.
#+plotThresh, results = "asis", fig.width=8, fig.height=6,echo=FALSE, message=FALSE, warning=FALSE
# Label points as fixations or saccades using the optimal velocity threshold
input <- input %>% mutate(fixation = ifelse(AV < optParam$value, 1, -1)) 
# Temporarily remove NAs to classify clusters
naIdx <- which(is.na(input$AV))
input$fixation[naIdx] = 0
input <- input %>% mutate(cluster = ifelse(fixation==1, cumsum(1-fixation)+1,0))
input$fixation[naIdx] = NA
# Arrange the clusters into ascending values
input <- input %>% group_by(trialIdx) %>%
  mutate(cluster = plyr::mapvalues(cluster, unique(sort(cluster)), 0:(length(unique(cluster))-1))) %>% ungroup()
# Plot the thresholded timecourse
start = round(min(na.omit(input$frameTime)), digits = 0)
end = round(max(na.omit(input$frameTime)), digits = 0)
timeSeq = breaks = seq(from = start, to = end, by = end/10)

# Filter the data to remove fixations shorter than 100ms
#clusterL <- input %>% filter(cluster>0) %>% group_by(trialIdx, cluster) %>% summarise(rleL = n()) %>% ungroup()
#shortCL <- which(clusterL$rleL<((Hz/2)/10))
clusterL <- input %>% filter(cluster>0) %>% group_by(trialIdx, cluster) %>% summarise(duration = max(frameTime) - min(frameTime)) %>% ungroup()
shortCL <- which(clusterL$duration<0.1)

shortIdx <- which(interaction(input$trialIdx, input$cluster) %in% interaction(clusterL$trialIdx[shortCL], clusterL$cluster[shortCL]))
if (length(shortIdx > 0)){
  input$cluster[shortIdx] = NA
  clusterSeq = seq(from = start, to = end, by = end/(end*2))
  # Code fixations <100ms as 0 
  input$fixation[shortIdx] = 0
}

for (i in 1:length(nTrials)) {
  current <- nTrials[i]
  cat(sprintf(subheader, runIdx, current))
  subdf <- input %>% filter(trialIdx==current)
  if (length(unique(na.omit(subdf$fixation)))==1){
    colScale = "#4292C6"
  } else (colScale = c("#4292C6", "#08306B", "#6BAED6"))
  
  xPlot <- ggplot(subdf) + 
    geom_line(aes(x = frameTime, y = Xf, 
                  colour = fixation, alpha = fixation), size = 1) +
    labs(x = "", y = "X-axis") +
    scale_colour_gradientn(colours = colScale, values = seq(-1, 1, by = 1)) +
    scale_x_continuous(breaks = timeSeq) +
    scale_alpha_continuous(range = c(1,1), na.value = 0) +
    theme_bw() + theme(legend.position = "none", 
                       panel.grid.minor.x = element_blank())
  yPlot <- ggplot(subdf) + 
    geom_line(aes(x = frameTime, y = Yf, colour = fixation, alpha = fixation), size = 1) +
    labs(x = "", y = "Y-axis") +
    scale_colour_gradientn(colours = colScale, values = seq(-1, 1, by = 1)) +
    scale_x_continuous(breaks = timeSeq) +
    scale_alpha_continuous(range = c(1,1), na.value = 0) +
    theme_bw() + theme(legend.position = "none")
  avPlot <- ggplot(subdf) + 
    geom_line(aes(x=frameTime, y = AV)) +
    labs(x = "Time elapsed", y = "Angular Velocity") + 
    scale_x_continuous(breaks = timeSeq) +
    geom_hline(yintercept = optParam$value, color = "red", linetype = "dashed") +
    theme_bw() + theme(legend.position = "none")
  # AVplot currently draws a line through points with no data present, might fix later
  xyvPlot <- plot_grid(xPlot, yPlot, avPlot, ncol = 1, align = "v")
  plot(xyvPlot)
  cat('\n\n')
}
#'
#' # Eyegaze Classification
#' Plots show saccades and fixations for each trial. Any fixations shorter than 100ms were discarded. The size of each point is
#' representative of the angular velocity and colour indicates its point in time.
#'  The bounding box represents the size and position of the image that was shown.
#' 
#+plotCluster, results = "asis", fig.width=8, fig.height=6,echo=FALSE, message=FALSE, warning=FALSE

# Filter the data to remove fixations shorter than 100ms
# clusterL <- input %>% filter(cluster>0) %>% group_by(trialIdx, cluster) %>% summarise(rleL = n()) %>% ungroup()
# shortCL <- which(clusterL$rleL<((Hz/2)/10))
# shortIdx <- which(interaction(input$trialIdx, input$cluster) %in% interaction(clusterL$trialIdx[shortCL], clusterL$cluster[shortCL]))
# input$cluster[shortIdx] = NA
#clusterSeq = seq(from = start, to = end, by = end/(end*2))

# size of the stimuli was 768*768 pixels and was fixed when presented
imgGrid <- data.frame(x = c((res_x/2)-(768/2), (res_x/2)+(768/2), (res_x/2)+(768/2), (res_x/2)-(768/2), (res_x/2)-(768/2)),
                y = c((res_y/2)-(768/2), (res_y/2)-(768/2), (res_y/2)+(768/2),(res_y/2)+(768/2), (res_y/2)-(768/2)))

# Identify timepoints that are outside of image ROI and calculate percentage spent gazing outside of ROI
input <- input %>% mutate(inROI = ifelse(Xf < min(imgGrid$x) | Xf > max(imgGrid$x) |
                                           Yf < min(imgGrid$y) | Yf > max(imgGrid$y), 
                                         0, 1)) %>%
  group_by(trialIdx) %>% mutate(roiValid = sum(inROI)/sum(!is.na(inROI)))


# Plot the final fixations
for (i in 1:length(nTrials)) {
  current <- nTrials[i]
  cat(sprintf(subheader, runIdx, current))
  subdf <- input %>% filter(trialIdx==current & !is.na(fixation) & fixation!=0)
  if (nrow(subdf)>1){
  xyFSplot <- ggplot(subdf) +
    geom_point(aes(x = Xf, y = Yf, size = AV, colour=frameTime), alpha = 0.5) +
    geom_path(aes(x = Xf, y = Yf),alpha = 0.1) +
    geom_polygon(data=imgGrid, aes(x=x, y=y), 
                 colour = "black", fill=NA, linetype="dashed", alpha = 0.05) +
    scale_x_continuous(limits = c(1,res_x), expand = c(0,0)) +
    scale_y_continuous(limits = c(1, res_y), expand = c(0,0)) +
    facet_grid(fixation ~ .,
               labeller = labeller(fixation = c("-1" = "saccades", "1" = "fixations"))) +
    labs(x = "X-Axis", y = "Y-axis") + theme_bw() + guides(size = FALSE) +
    scale_colour_gradientn(name = "time elapsed", breaks = clusterSeq, labels = as.character(clusterSeq) ,colours=rainbow(4)) +
    theme(legend.key.height = unit(0.5,'inch'))
  plot(xyFSplot)}
  cat('\n\n')
}
