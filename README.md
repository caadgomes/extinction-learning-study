# extinction-learning-study

A collection of scripts used in the manuscript "Predicting individual differences of fear and cognitive learning and extinction"

########## GENERAL INFORMATION #############

All scripts are provided in the folder scripts. Each script has an associated .txt file with a summary of its functionality, instructions how to run it, and which libraries and files are required.

############ ORDER OF SCRIPTS ###############

The scripts should be run in the following order:

01. pipeline_BIDS.py
This script will convert all files (functional, structural, diffusion, physiological, behavioural) into an appropriate BIDs format.

02. pipeline_fmriprep.py
This script will run fmriprep on BIDs functional and anatomical data, and create the corresponding preprocessed files.

03. pipeline_denoising.py
This script performs denoising of the preprocessed functional images from step 2.

04. pipeline_ROIs.py
This script uses the preprocessed anatomical images as well as corresponding freesurfer outputs from step 2 to create subject-specific ROIs.

05. pipeline_suit.m
This script uses the preprocessed anatomical images as well as corresponding freesurfer outputs from step 2 to create subject-specific cerebellum ROIs.

06. pipeline_FC.py
This script computes different types of functional connectivity metrics within selected ROIs (step 4) from the denoised functional images from step 3.

07. pipeline_SC.py
This script computes structural connectivity within selected ROIs (step 4) from the dwi images.

08. pipeline_EC.m
This script computes effective connectivity within selected ROIs (step 4) from the denoised functional images from step 3.

09. pipeline_EDA.m, pipeline_EDA_DCM.m, pipeline_EDA_SF
These scripts computes EDA learning metrics (based on DCM or SF) from EDA raw data stored in BIDs.

10. pipeline_learning.R
This script computes learning metrics either from EDA (step 9) or behavioural data within selected ROIs (step 4) from the denoised functional images from step 3.

11. pipeline_modelling.R
This script predicts learning (step 10) from brain connectivity measures (steps 6, 7 and 8).

12. pipeline_generalisability1.R
This script computes generalisation of predictions from step 11.

13. pipeline_LOGO.py
This script computes generalisation of predictions from step 11 using LOGO.

############## DEPENDENCIES #################

All scripts either require python, R or matlab to be installed (check specific scripts for details).

The python scripts (.py) require the latest versions of the following libraries:
ast, bctpy, collections, csv, datetime, functools, errno, gzip, itertools, json, multiprocessing, nibabel, nilearn, nipype, numpy, operator, os, pandas as pd, pathlib, pathlib2, psutil, pydicom, pycwt, re, scipy, shlex, shutil, sklearn, string, subprocess, tslearn, warnings (these libraries can be installed using: pip install [lib])

The R scripts (.r) require the latest version of the following libraries:
lme4, emmeans, ggplot2, ggcorrplot, dplyr, stringr, tibble, tidyr, lmerTest, car, caret, glmnet, doMC (these libraries can be installed using: install.libraries("lib")

The matlab scripts (.m) require SUIT, PsPM and SPM toolboxes.
