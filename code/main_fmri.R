#STEP 0: FOLDER STRUCTURE
#Make sure to have a clear organisation of files in place. 

#Workdir
#raw_files
#pp_1
#pp_2
#preprocessed_files
#pp_1
#pp_2
#html_reports
#pp_1
#pp_2

#STEP 1: LOAD LIBRARIES AND SET PATH
library(dplyr)
library(rmarkdown)
# Specify the paths to the data and the code
dataDir <- choose.dir() #parent folder to rawdata
workDir <- choose.dir() #dir containing the R code
rawdata <- paste0(dataDir, '/rawdata/')
#STEP 2: SPECIFY EXPERIMENT PARAMETERS
#Make sure to specify all the parameters necessary!

preprocDir <- paste0(dataDir, '/preproc/')
dir.create(preprocDir)
# Screen resolution in pixels as Eyelink returns pixel coordinates
res_x = 1920
res_y = 1080
scr_x = 368
scr_y = 202

Hz = 1000 # Sampling rate
win = 5 # Window length for median filter
distance = 67 # Viewing distance (in cm) which was fixed

#STEP 3: PREPROCESSING
subdirs <-  list.dirs(rawdata, recursive = F)

for (s in 1:length(subdirs)){
  subject = strsplit(subdirs,split = '/')[[s]][4]
  rawfiles <- list.files(path=paste0(rawdata, subject), pattern = c("_run|eyedata.csv)"))
  nFiles <- length(rawfiles)
  
  #trial types are (1,3) for perception of analog, digital
  #trial types are (2,4) for imagery of analog, digital
  #trial type 5 is the first trial per block
  #trial type 0 is fixation cross
  trialcode = c(0,1,2,3,4)
  
  #PREPROCESSING
  for (f in 1:nFiles){
    rawDF <- read.csv(paste0(rawdata, '/', subject, '/', rawfiles[f]), 
                      stringsAsFactors = FALSE, header = TRUE)
    allRuns = unique(rawDF$runIdx)
    pp = unique(rawDF$ppID)
    
    if (!dir.exists(file.path(preprocDir, pp))){
      dir.create(paste0(preprocDir,pp))
    }
    
    for (thisTrial in 1:length(trialcode)){
    #for (thisTrial in 1){
      outputFile = paste0(pp, "-run", allRuns, "-trialtype", trialcode[thisTrial], "-preproc.html")
      outDir = paste0(preprocDir, pp, "/report/")
      input <- rawDF %>% filter(runIdx==allRuns & trialCode==trialcode[thisTrial])
      
      #RUN THE DAMN THING
      rmarkdown::render(input = paste0(workDir, "/preproc_eyelink.R"),
                        output_file = outputFile, output_dir = outDir,
                        clean = TRUE)
      output = paste0(preprocDir, pp, "/csv/", '/', pp, "-run", allRuns, "-stim", trialcode[thisTrial],"-preproc.csv")
      dir.create(paste0(preprocDir, pp, "/csv/"))
      write.csv(input, output)
    }
  }
}

# GAZEPATH MAPPING

gazepathDir <- paste0(dataDir, '/gazepath/')
dir.create(gazepathDir)

# Set location of fMRI experiment playlists
pathtoplaylists <- choose.dir()

subdirs = list.dirs(preprocDir, recursive = F)

# Loop through data per subject
for (s in 1:length(subdirs)){
  subject = strsplit(subdirs,split = '/')[[s]][4]
  for (run in 1:4){
    preproc = list.files(path=paste0(preprocDir, subject,'/csv/'), pattern=paste0(subject, '-run', run))[-1] # the [-1] removes fixation cross output from the list
    preproc <- paste0(paste0(preprocDir,subject, '/csv/'), preproc)
    if (length(preproc)>0){
      for (t in c(1,3)){
        # Only run if both have been processed
        input1 <- read.csv(preproc[t], stringsAsFactors = F, header = T)[-1]
        input2 <- read.csv(preproc[t+1], stringsAsFactors = F, header = T)[-1]
        outDir = paste0(gazepathDir, subject, "/report/")
        outputFile = paste0(subject, "_stim",t, '_run', run, "-gazepath.html")
        # Read playlist
        playlists  <- list.files(path=pathtoplaylists, pattern='pilot_example')
        playlists <- paste0(pathtoplaylists,'/', playlists)
        playlist <- read.csv(playlists[run])
        if (exists('report')){remove(report)}
        rmarkdown::render(input = paste0(workDir, "/gazepath_eyelink.R"),
                          output_file = outputFile, output_dir = outDir,
                          clean = TRUE)
        output = paste0(gazepathDir, subject,'/csv/', subject, "-run_", run, "-stim_", t ,"-gazepath.csv")
        dir.create(paste0(gazepathDir, subject, '/csv/'))
        write.csv(report, output)
      }
    }
  }
}

#MATRIX SIMILARITY

similarityDir <- paste0(dataDir, "/similarity/")
dir.create(similarityDir)

# CALCULATED ACROSS RUNS RATHER THAN FOR EACH RUN
for (s in 1:length(subdirs)){
  subject = strsplit(subdirs,split = '/')[[s]][4]
  subdata = paste0(preprocDir, subject, '/csv/')
  fxFiles <- list.files(path = subdata, pattern = 'stim0')
  fxFiles <- paste0(subdata, fxFiles)
  concFX <- do.call(rbind, lapply(fxFiles, read.csv))
  
  for (i in c(1,3)){
    stimFiles <- list.files(path = subdata, pattern = paste0('stim', i))
    stimFiles <- paste0(subdata, stimFiles)
    concStim = do.call(rbind, lapply(stimFiles, read.csv))
    ImgFiles <- list.files(path = subdata, pattern = paste0('stim', i+1))
    ImgFiles <- paste0(subdata, ImgFiles)
    concImg = do.call(rbind, lapply(ImgFiles, read.csv))

    # Read playlist
    playlists  <- list.files(path=pathtoplaylists, pattern='pilot_example')
    playlists <- paste0(pathtoplaylists,'/', playlists)
    playlist <- do.call(rbind, lapply(playlists, read.csv))
    
    stims <- playlist$trialNo[playlist$stimFile %in% playlist$imgClock & playlist$trialCode==i]
    # run each combo 0-i 0-i+1 i-i+1
    for (j in 0:2){
      if (j==0){
        combo <- paste0(0,i)
        input1 <- concFX
        input2 <- concStim
      } else if (j==1){
        combo <- paste0(0,i+1)
        input1 <- concFX
        input2 <- concImg
      } else if (j==2){
        combo <- paste0(i,i+1)
        input1 <- concStim[concStim$trialIdx %in% stims,]
        input2 <- concImg
      }
      if (i==1){
        input1$clockCode=1
        input2$clockCode=1
      } else {
        input1$clockCode=2
        input2$clockCode=2
      }
      outDir = paste0(similarityDir, subject, "/report/")
      outputFile = paste0(subject, "-allRuns-stimpair_", combo ,"-Rv_corr.html")
      rmarkdown::render(input = paste0(workDir, "/RV_eyelink.R"),
                        output_file = outputFile, output_dir = outDir,
                        clean = TRUE)
      if (!exists('rvDF')){
        rvDF <- as.data.frame(rvOut)
      } else {
        rvDF <- rbind(rvDF, rvOut)
      }
    }
  }
}
dir.create(paste0(similarityDir, subject, "/csv/"))
write.csv(rvDF, paste0(similarityDir, subject, '/csv/Rv_all_matched.csv'))



## Similarity within each run ####
# for (s in 1:length(subdirs)){
#   subject = strsplit(subdirs,split = '/')[[s]][7]
#   subdata = paste0(preprocPath, '/' ,subject, '/')
#   for (r in 1:4){
#     thisRun = paste0("run", r)
#     pat = paste0(subject, "-", thisRun)
#     inputFiles <- list.files(path = subdata, pattern = pat)
#     if (length(inputFiles)>0){
#       for (i in c(1,3)){
#         input1 <- read.csv(paste0(subdata, inputFiles[i]), stringsAsFactors = F, header = T)
#         input2 <- read.csv(paste0(subdata, inputFiles[i+1]), stringsAsFactors = F, header = T)
#         if (i==1){
#           input1$clockCode=1
#           input2$clockCode=1
#         } else {
#           input1$clockCode=2
#           input2$clockCode=2
#         }
#         outDir = paste0(dataDir, "reports/", subject)
#         outputFile = paste0(subject, "_", thisRun, "_stim", i ,"-Rv_corr.html")
#         rmarkdown::render(input = paste0(workDir, "code/fixation_similarity.R"),
#                           output_file = outputFile, output_dir = outDir,
#                           clean = TRUE)
#         #bind the output to a dataframe
#         if (!exists('rvDF')){
#           rvDF <- as.data.frame(rvOut)
#         } else {
#           rvDF <- rbind(rvDF, rvOut)
#         }
#       }
#     }
#   }
# }