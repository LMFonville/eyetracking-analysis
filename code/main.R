# STEP 0: FOLDER STRUCTURE
# Make sure to have a clear organisation of files in place. 
# I use separate folders for the raw data, the preprocessed data, and the html reports

#Workdir
  #raw_files
  #preprocessed_files
    #pp_1
    #pp_2
  #html_reports
    #pp_1
    #pp_2

#STEP 1: LOAD DATA & LIBRARIES
library(dplyr)
library(rmarkdown)
# Let's load the datafile which contains all the trials
if (!dir.exists("C:/Users/fonvill/")){
  workDir <-  "C:/Users/leonf/Dropbox/1.dataviz/eyetracking-analysis/"
} else{
  workDir <-  "E:/1.Work/visual-imagery-clocktask/"
  
}
rawDF <- read.csv(paste0(workDir,"rawdata/", "exampleData.csv"), 
                  stringsAsFactors = FALSE, header = TRUE)
pp = unique(rawDF$ppID)

# If directory to store preprocessed data does not exist make the output directory
preprocDir = paste0(workDir, "preproc/pp_", pp, "/")
if (!dir.exists(preprocDir)){
  dir.create(preprocDir)
}
# The directory to store reports will be created automically when rendering if it does not exist
reportDir = paste0(workDir, "reports/pp_", pp) 

# STEP 2: SPECIFY YOUR EXPERIMENT PARAMETERS. (NON-OPTIONAL)
# Make sure to specify all the parameters necessary!
# Hz is the sampling rate at which you acquired your data
Hz = 120
# Screen dimensions in mm
scr_x = 509
scr_y = 286

# MY DATA IS ACQUIRED IN THE ACTIVE COORDINATE SYSTEM SO EYEGAZE DATA RANGES FROM 0-1. 
# IN ORDER TO CALCULATE THE ANGLE I NEED TO GET THE DISTANCE BETWEEN POINTS IN MM.
# IF YOUR DATA IS ALREADY CONVERTED TO DISTANCE IN MM CHANGE THESE VALUES TO 1. 

# STEP 3: FILTER YOUR DATASET (OPTIONAL)
# If you only want to use a subset of your data, certain types of trials or conditions, then use filtering to create this smaller dataframe
# I've got multiple runs and trial types so I am going to loop through my runs (1-6) and different trial types (1 & 3) these to produce reports
allRuns = unique(rawDF$runIdx)
trials = c(1,3)

for (thisRun in 1:length(allRuns)){
  for (thisTrial in 1:length(trials)){
    outputFile = paste0("pp",pp, "-run", allRuns[thisRun], "-stim", trials[thisTrial], "-preproc.html")
    input <- rawDF %>% filter(runIdx==allRuns[thisRun] & trialCode==trials[thisTrial])
    # STEP 4: RUN THE DAMN THING
    rmarkdown::render(input = paste0(workDir, "code/preproc.R"),
                      output_file = outputFile, output_dir = reportDir,
                      clean = TRUE)
    output = paste0(preprocDir,"pp", pp, "-run", allRuns[thisRun], 
                    "-stim", trials[thisTrial],"-preproc.csv")
    write.csv(input, output)
  }
}

#STEP 5: CHECK FOR THE FILES IN YOUR SPECIFIED FOLDER
