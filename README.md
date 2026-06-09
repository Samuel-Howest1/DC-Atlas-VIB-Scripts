# DC-Atlas-VIB-Scripts

This is a repository for my interschip in the VIB creating a DC cell Atlas. This repository contains all scripts related to this project.

This project use renv to keep package versions consitents between devices

# Content

The 3 main folders (Reports,Rshiny tool,Scripts) contain the following:

## Scripts

the Scripts folder contains 2 main sub-folders: "Processing of Data" and "Integration and Multimodel".

The "Processing of Data" folder contains all script used for the reports. It contains 3 scripts in total: a Pre-processing script, that filters out outlieres cells, finds doublets and annotates cells. a Hashing script that identifies each hashtag and find the doublets and negatives. and finally a Post-processing script that is used to manually correct unknown/suspicious cells, either correcting the annotation or removing doublets.

The "Integration and Multimodel" folder script to merge all resulting data of the post-processing script, scripts for all 3 integration methods used, a script to score these integration's ,A script to create feature plots of biomarks for the CDC1 cell life cycle, and finally a script to make the data weighted multi modal.

## Reports

The Reports folder contains all the reports for the processed datasets. This means only what was done to processes the data. The datasets and the biological result are not present. Additionally not all are published and are kept private by the VIB until the paper is complete.


## Rshiny tool
