## Introduction

This repository describes a series of processing steps and some analyses that I've used on some eye-tracking data acquired during fMRI. These can be run directly from R and will produce a csv output of the preprocessed data and will also render reports as html files for quality control. There are several R packages needed to run this pipeline (currently: dplyr, ggplot2, cowplot, markdown). Below are links of examples of the output as an html report with figures and code. 

This data was quite messy and required downsampling, smoothing, and clipping to become useable. The main motivation for this repository was to make the processing steps easier to follow by visualising what is happening and providing a direct reference to the underlying code in the report. This project is no longer ongoing as I don't have time to pursue it further.

## Report examples

These are example runthroughs from a pilot scan. The html previews don't have the same functionality as the reports (code folding, table of contents, layout).

Use CTRL-click to open in a new tab! 

Preprocessing
* [Example 1](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/preproc/pilot/report/pilot-run1-trialtype1-preproc.html)

* [Example 2](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/preproc/pilot/report/pilot-run1-trialtype2-preproc.html)

Gazepath
* [Example 1](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/gazepath/pilot/report/pilot_stim1_run1-gazepath.html)

* [Example 2](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/gazepath/pilot/report/pilot_stim3_run1-gazepath.html)

Similarity
* [Example 1](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/similarity/pilot/report/pilot-allRuns-stimpair_12-Rv_corr.html)

* [Example 2](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/similarity/pilot/report/pilot-allRuns-stimpair_01-Rv_corr.html)
