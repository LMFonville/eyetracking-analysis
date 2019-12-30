#' ---
#' title: Gazepath Similarity
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
#' 
#' This report examines similarity in the sequences of eyegaze fixations between
#' two trial types. This is mostly based on the work from 
#' [Mathot et al. (2012)](https://bop.unibe.ch/JEMR/article/view/2326).
#' For convenience we will refer to the gazepath for the stimulus trials as X
#' and for the imagery trials as Y.
#' 
#+LoadLib, echo=FALSE, warning=FALSE, message=FALSE
library(dplyr)
library(ggplot2)
library(cowplot)

#' # Summary
#+prepSum
#Some housekeeping
ppID = unique(c(input1$ppID,input2$ppID))
runIdx = unique(c(input1$runIdx, input2$runIdx))
stimcode <- unique(input1$trialCode)
imgcode <- unique(input2$trialCode)

if (stimcode==1 & imgcode==2){stimType = "Analog Clocks"
} else if(stimcode==3 & imgcode==4){stimType = "Digital Clocks"
} else {stimType = "Invalid Code"}

#prep for reporting trials
subheader <- "### Run %d Stim Trial %d Img Trial %d

"

# Get trial index for stim and img trials
stims <- playlist$trialNo[playlist$stimFile %in% playlist$imgClock & playlist$trialCode==stimcode]
imgs <- playlist$trialNo[playlist$imgClock %in% playlist$stimFile & playlist$trialCode==imgcode]
imglist <- list()
for (stimL in 1:length(stims)){
  imglist[[stimL]] <- playlist$trialNo[playlist$imgClock %in% playlist$stimFile[stims[stimL]]]
}
# Remove stim trials with no matching img trial
input1 <- input1[input1$trialIdx %in% stims,]
input2 <- input2[input2$trialIdx %in% imgs,]


#'  |                   |              |
#'  |:------------------|-------------:|
#'  | Participant ID:   | `r ppID`     |
#'  | Experiment Run:   | `r runIdx`   | 
#'  | Stimulus Type:    | `r stimType` | 

#' # Calculations
#+ calcDistance, warning=FALSE
#' The mean of each fixation cluster was used as cluster centroid and timestamps were scaled to improve distance calculations. 
#' The Euclidean distance was then calculated for every centroid mapped to its 
#' nearest point in the other scanpath and the distance is normalised by
#' dividing the sum of Euclidean distances by the number of fixations in
#' whichever is the longest sequence. Here, we calculate a 3-dimensional
#' Euclidean distance using the XY coordinates as well as the timepoint.
#' A greater similarity is represented by a lower distance.

input1C <- input1 %>% filter(cluster>0) %>% 
  group_by(trialIdx, cluster) %>% 
  summarise(Xc = mean(Xf), Yc = mean(Yf), Fc = mean(frameTime)) %>% 
  group_by(trialIdx) %>% mutate(zTime = (Fc - min(Fc))/(max(Fc)-min(Fc)))
input2C <- input2 %>% filter(cluster>0) %>% 
  group_by(trialIdx, cluster) %>% 
  summarise(Xc = mean(Xf), Yc = mean(Yf), Fc = mean(frameTime)) %>% 
  group_by(trialIdx) %>% mutate(zTime = (Fc - min(Fc))/(max(Fc)-min(Fc)))
input1C$zTime <- ifelse(is.na(input1C$zTime), 0.5, input1C$zTime)
input2C$zTime <- ifelse(is.na(input2C$zTime), 0.5, input2C$zTime)
gazepath.Fun <- function(x, y, stimcode){
  #calculate distances
  Dm = matrix(data=NA, nrow = length(x$Xc), ncol = length(y$Xc))
  for (i in 1:length(x$Xc)){
    for (j in 1:length(y$Xc)){
      # 3 dimensional Euclidean Distance
      Dm[i,j] = sqrt((x$Xc[i] - y$Xc[j])^2 + (x$Yc[i] - y$Yc[j])^2 + (x$zTime[i] - y$zTime[j])^2)
    }
  }
  xMin <- apply(Dm, 1, function(x) min(x)) #lowest D per row
  xMinIdx <- apply(Dm, 1, function(x) which.min(x)) #index of lowest D per row
  yMin <- apply(Dm, 2, function(x) min(x)) #lowest D per column
  yMinIdx <- apply(Dm, 2, function(x) which.min(x)) #index of lowest D per column
  #Add index to xy dataframe
  xy <- rbind(x,y)
  xy$trialCode <- rep(c(stimcode, stimcode+1), c(nrow(x), nrow(y)))
  xy$Didx <- c(yC[xMinIdx], xC[yMinIdx])
  #Calculate the total distance
  Dst <- sum(xMin,yMin)/max(length(xMin), length(yMin)) #Total sum of distances given number of clusters
  
  # Make a df for plotting
  xy <- xy %>% mutate(label = ifelse(trialCode==stimcode, 1, 2), group = 1:n())
  plotdf <- xy
  limit <- nrow(plotdf)
  for (c in 1:limit){
    cluster = xy$Didx[c]
    if (xy$trialCode[c]==stimcode){
      #find the matching Y cluster for the X centroid
      rowIdx = which(xy$cluster==cluster & xy$trialCode==stimcode+1)
    } else if (xy$trialCode[c]==stimcode+1){
      #find the matching X cluster for the Y centroid
      rowIdx = which(xy$cluster==cluster & xy$trialCode==stimcode)
    }
    plotdf <- rbind(plotdf, xy[rowIdx,])
  }
  #add variable to indicate which rows are new and adjust the label for new rows
  plotdf$new <- c(rep(0, limit), rep(1, limit))
  plotdf$label[plotdf$new==1] <- ifelse(plotdf$label[plotdf$new==1]==1, 2, 1)
  #replace group values so that new rows match their mapped centroids
  clusters <- sort(unique(plotdf$Didx))
  for (i in 1:length(clusters)){
    centroid = clusters[i]
    #update the group value for each added row to match the new centroid with the original cluster
    oldG <- which(plotdf$cluster==centroid & plotdf$new==1)
    newG <- which(plotdf$Didx==centroid & plotdf$new==0)
    plotdf$group[oldG] <- plotdf$group[newG]
  }
  #rename label for plotting and add a shape variable to improve differentiating between centroids
  plotdf$label <- factor(plotdf$label, levels = c(1,2), labels = c("X mapped to Y", "Y mapped to X"))
  plotdf$plotShape <- interaction(plotdf$trialCode, plotdf$new)
  
  
  # Make trialreport and add to a dataframe
  trialreport <- data.frame(stimtrial = unique(x$trialIdx),
                            imgtrial = unique(y$trialIdx),
                            D = Dst,
                            Xpoints = length(unique(x$Xc)),
                            Ypoints = length(unique(y$Xc)),
                            XtoY = length(unique(yC[xMinIdx])),
                            YtoX = length(unique(xC[yMinIdx])))
  
  output <- list()
  output$report <- trialreport
  output$plotdf <- plotdf
  return(output)
}


#' # Point-to-Point Mapping
#' Plots show the original scanpath, the nearest centroid each fixation cluster is
#' mapped to and the scanpath using those mappings.
#+PlotPaths, results="asis", fig.height=10, fig.width=12, echo=FALSE, warning=FALSE, message=FALSE
for (stimtrial in 1:length(stims)){
  x <- input1C %>% filter(trialIdx==stims[stimtrial])
  xC <- input1C %>% filter(trialIdx==stims[stimtrial]) %>% pull(cluster)
  for (imgtrial in 1:length(imglist[[stimtrial]])){
    img_idx <- imglist[[stimtrial]][imgtrial]
    y <- input2C %>% filter(trialIdx==img_idx)
    yC <- input2C %>% filter(trialIdx==img_idx) %>% pull(cluster)
    if (length(xC)>0 & length(yC)>0){
      gazepath <- gazepath.Fun(x, y, stimcode)
      if (!exists('report')){report <- gazepath$report} else {
        report <- rbind(report, gazepath$report)}
      start = round(min(gazepath$plotdf$zTime), digits = 0)
      end = round(max(gazepath$plotdf$zTime), digits = 0)
      fixSeq = seq(from = start, to = end, by = end/(end*2))
      cat(sprintf(subheader, runIdx, stims[stimtrial], img_idx))
      pathP <- ggplot(gazepath$plotdf %>% filter(new==0), aes(x=Xc, y = Yc)) + 
        geom_point(aes(shape = factor(trialCode), 
                       colour = zTime), size = 3) + 
        geom_path() + 
        facet_wrap(~ label, labeller = 
                     labeller(label = c("X mapped to Y" = "X gazepath", 
                                        "Y mapped to X" = "Y gazepath" ))) + 
        labs(x = "X", y = "Y", colour = "time") +
        scale_x_continuous(limits = c(0, res_x), expand=c(0,0)) + 
        scale_y_continuous(limits = c(0, res_y), expand=c(0,0)) +
        theme_bw() + 
        scale_colour_gradientn(name = "time elapsed", breaks = fixSeq, labels = as.character(fixSeq) ,colours=rainbow(4)) +
        scale_shape_manual(values = c(15, 19), guide = FALSE) +
        theme(legend.position = "none")
      mapP <- ggplot(gazepath$plotdf, aes(x=Xc, y=Yc)) + 
        geom_point(aes(shape = factor(plotShape), 
                       colour = zTime), size = 3) +
        geom_line(aes(group = group)) + facet_wrap(~ label) +
        labs(x = "X", y = "Y", colour = "time") +
        scale_x_continuous(limits = c(0, res_x), expand=c(0,0)) + 
        scale_y_continuous(limits = c(0, res_y), expand=c(0,0)) +
        theme_bw() + 
        scale_colour_gradientn(name = "time elapsed", breaks = fixSeq, labels = as.character(fixSeq) ,colours=rainbow(4)) +
        scale_shape_manual(values = c(22, 21, 15, 19)) +
        theme(legend.position = "none")
      newP <- ggplot(gazepath$plotdf, aes(x=Xc, y=Yc)) + 
        geom_point(aes(shape = factor(trialCode), 
                       colour = zTime,
                       alpha=new), size = 3) +
        geom_path(data=gazepath$plotdf %>% filter(new==1)) + 
        facet_wrap(~ label, labeller = 
                     labeller(label = c("X mapped to Y" = "New X-mapped gazepath",
                                        "Y mapped to X" = "New Y-mapped gazepath" ))) +
        labs(x = "X", y = "Y", colour = "time elapsed (seconds)") +
        scale_x_continuous(limits = c(0, res_x), expand=c(0,0)) + 
        scale_y_continuous(limits = c(0, res_y), expand=c(0,0)) +
        theme_bw() + 
        scale_alpha_continuous(range = c(0,1), na.value = 0, guide = FALSE) +
        scale_shape_manual(values = c(19, 15), guide = FALSE) + 
        scale_colour_gradientn(name = "time elapsed", breaks = fixSeq, labels = as.character(fixSeq) ,colours=rainbow(4)) +
        theme(legend.position = "bottom")
      xyP <- plot_grid(pathP, mapP, newP, ncol = 1, nrow = 3)
      plot(xyP)
      cat('\n\n')
    }
  }
}

#' # Distance Summary
#' 
knitr::kable(report, digits = 2, col.names = c("Stimulus Trial", "Imagery Trial", "Distance", "points X", "points Y", "X to Y", "Y to X"))
