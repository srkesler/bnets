#!/bin/bash

# Shelli Kesler 6/17/15
#This script copies the ROI folder into the subject directory, applies the transformation 
#to each ROI, removes the original ROIs, creates a path text file for the ROIs, computes
#number of tracks and track stats for each pair of ROIs and the extracts the volume for
#each ROI. 

export PATH=$PATH:/path/umcp_13_dev #edit: path to UMCP folder
# http://ccn.ucla.edu/wiki/index.php/UCLA_Multimodal_Connectivity_Package

SUBJS=$(cat /Volumes/RAID24/SJ_DTI/paths.txt) #edit path to paths file
shopt -s extglob #turn on extended pattern matching to use negative wildcards (!)

#copy ROI folder into subject directory
for subjdir in ${SUBJS}
do 
	cp -r /path/ROIs ${subjdir} #path to ROI folder

#apply transformation to each ROI to create DTI native space ROI files

	cd ROIs
	rois=`ls *.gz`
	pre='w'

		for i in ${rois}
		do
    		flirt -in ${i} -applyxfm -init ../FSPGR_brain_B0_MNI_inv.mat -out ${pre}${i} \
    		-paddingsize 0.0 -interp trilinear -ref ../nodif_brain.nii.gz
		done
	
	rm !(w*) #delete all files that do not begin with w
	
#create ROIpath.txt for each subject

	rois=`ls -d -1 $PWD/*.gz` 
		for f in ${rois}
		do
			echo "$f" >> ../ROIpaths.txt 
		done

#count number of tracks between ROIs and calculate track FA stats

	cd ../
	run_tracks.py -t dti.trk -m ROIpaths.txt -o subj -c -s --statimg=dti_fa.nii 

#extract volumes in mm3 for ROIs
	cd ROIs
	rois2=`ls -d -1 *.gz` 
		for reg in ${rois2}
		do
			fslstats ../FSPGR_brain_B0.nii.gz -k ${reg} -V | awk '{print $2}' >> ../ROIvols.csv
		done
	
done

shopt -u extglob #turn off extended pattern matching
