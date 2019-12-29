## Introduction

This repository describes a series of processing steps that I've used to analyse eye-tracking data acquired during fMRI. These can be run directly from R and will produce a csv output of the preprocessed data and will also render reports as html files for quality control. There are several R packages needed to run this pipeline (currently: dplyr, ggplot2, cowplot, markdown). Below are links of examples of the output as an html report with figures and code. This data was quite messy and required downsampling, smoothing, and clipping to become useable. 

The main motivation for this repository was to make the processing steps easier to follow by visualising what is happening and providing a direct reference to the underlying code in the report. 

## Report examples

These are examples from my data that have problematic features to showcase how the pipeline deals with both valid and invalid data.

Use CTRL-click to open in a new tab! 

* [Example 1](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/reports/pp_999/pp999-run1-stim1-preproc.html)

* [Example 2](http://htmlpreview.github.io/?https://github.com/LMFonville/eyetracking-analysis/blob/master/reports/pp_999/pp999-run1-stim3-preproc.html)

