#!/bin/bash

#Shelli Kesler 5/7/14

#change path below to location of Diffusion Toolkit
export PATH=$PATH.:/Applications/Diffusion\ Toolkit.app/Contents/MacOS

#This script performs DTI tensor reconstruction and tractography and transforms ROIs 
#into DTI native space

# Check the dti.trk file and define ROIs before moving on to dcon3

SUBJS=$(cat /Volumes/EHD/Study/paths.txt) #edit path to paths file

for subjdir in ${SUBJS}
do 
	cd ${subjdir}

	# Mask the eddy corrected DTI data with the B0 mask
	fslmaths dti_corr.nii.gz -mas nodif_brain_mask.nii.gz dti_corr_masked
	
	#Tensor reconstruction and deterministic tractography
	#Ensure gradient table txt file path and bvals are correct here
	dti_recon "dti_corr_masked.nii.gz" "dti" -gm "/path/BVECS.txt" \
		-b 850 -b0 auto -p 3 -sn 1 -ot nii 
	
	#change track angle here if necessary (default = 60 degrees)
	dti_tracker "dti" "track_tmp.trk" -at 60 -m "dti_dwi.nii" -it nii

	spline_filter "track_tmp.trk" 1 "dti.trk"
	
	#Extract the brain from the T1 weighted anatomic volume
	bet FSPGR.nii.gz FSPGR_brain
	
	#Co-register the T1 to the skull stripped B0 to create a T1 image that is in DTI space
	flirt -in FSPGR_brain.nii.gz -ref nodif_brain.nii.gz -out FSPGR_brain_B0.nii.gz \
	-omat FSPGR_brain_B0.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 \
	-searchrz -90 90 -dof 12  -interp trilinear

	#Co-register the DTI space T1 to MNI space
	flirt -in FSPGR_brain_B0.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz \
	-out FSPGR_brain_B0_MNI.nii.gz -omat FSPGR_brain_B0_MNI.mat -bins 256 -cost corratio -searchrx \
	-90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear
	
	#Save the inverse of the above transformation
	convert_xfm -omat FSPGR_brain_B0_MNI_inv.mat -inverse FSPGR_brain_B0_MNI.mat
	
	#unzip the DTI space T1 (Rex requires uncompressed version)
	#gunzip -k FSPGR_brain_B0.nii.gz

done
