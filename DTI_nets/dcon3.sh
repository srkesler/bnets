#!/bin/bash

# Shelli Kesler 6/17/15
#This script copies the ROI folder into the subject directory, applies the transformation 
#to each ROI, removes the original ROIs, creates a path text file for the ROIs, computes
#number of tracks and track stats for each pair of ROIs and the extracts the volume for
#each ROI. 

#This script must be inside the STUDY directory.

export PATH=$PATH:/path/umcp_13_dev #edit: path to UMCP folder
# http://ccn.ucla.edu/wiki/index.php/UCLA_Multimodal_Connectivity_Package

STUDY='/path/' #path to study directory containing all subject folders
SUBJS=`ls -d */`
datadir='connectome/' #directory containing preprocessed DTI and anatomic data
datadir2='connectome/ROIs/' #ROI directory within subject's folder after it is copied there
shopt -s extglob #turn on extended pattern matching to use negative wildcards (!)

#copy ROI folder into subject directory
for subjdir in ${SUBJS}
do 
	currentdir=${STUDY}${subjdir}${datadir}
	cp -r /path/ROIs ${currentdir} #path to ROI folder
done

#apply transformation to each ROI to create DTI native space ROI files
for subjdir in ${SUBJS}
do 
	currentdir=${STUDY}${subjdir}${datadir2}
	cd ${currentdir}
	rois=`ls *.gz`
	pre='w'

	for i in ${rois}
	do
    	flirt -in ${i} -applyxfm -init ../FSPGR_brain_B0_MNI_inv.mat -out ${pre}${i} \
    	-paddingsize 0.0 -interp trilinear -ref ../nodif_brain.nii.gz
	done
	
	rm !(w*) #delete all files that do not begin with w
done
	
#create ROIpath.txt for each subject
for subjdir in ${SUBJS}
do 
	currentdir=${STUDY}${subjdir}${datadir2}
	cd ${currentdir}
	rois=`ls -d -1 $PWD/*.gz` 
		for f in ${rois}
		do
			echo "$f" >> ../ROIpaths.txt 
		done
	
done

#count number of tracks between ROIs and calculate track FA stats
for subjdir in ${SUBJS}
do 
	currentdir=${STUDY}${subjdir}${datadir}
	cd ${currentdir}
	run_tracks.py -t dti.trk -m ROIpaths.txt -o subj -c -s --statimg=dti_fa.nii 
done

#extract volumes in mm3 for ROIs
for subjdir in ${SUBJS}
do 
	currentdir=${STUDY}${subjdir}${datadir2}
	cd ${currentdir}
	rois2=`ls -d -1 *.gz` 
		for reg in ${rois2}
		do
			fslstats ../FSPGR_brain_B0.nii.gz -k ${reg} -V | awk '{print $2}' >> ../ROIvols.csv
		done
	
done

shopt -u extglob #turn off extended pattern matching
