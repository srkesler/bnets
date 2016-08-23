#!/bin/bash

#Shelli Kesler 12/8/14 
# Check mask files before moving on to dcon2

SUBJS=$(cat /Volumes/EHD/Study/paths.txt) #edit path to paths file

for subjdir in ${SUBJS}
do 
	cd ${subjdir}
	
	#eddy current correction of raw DTI data
	eddy_correct data.nii.gz dti_corr 0

	#extract the first B0 volume
	fslroi dti_corr.nii.gz nodif 0 1
	
	#skull strip the B0 volume and create a binary mask of it
	bet nodif.nii.gz nodif_brain -f 0.5 -g 0 -m
done
