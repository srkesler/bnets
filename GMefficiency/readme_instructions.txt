1. Arrange your subject folders accordingly and have filelist.txt ready as shown in the folder.

2. Next run dicm2nii_conversion.m which will convert dicoms to nifti files and also will name the .nii file accordingly.

3. Next run dartelized vbm using spm8 gui. There is a sample dartelvbm_71_job.mat for example.

4. After that is ready make a filelist_GM.txt and paths.txt as shown in the folder.

5. Next run nodalefficiencies_calculation.m and that will save nodalefficiency_stats.mat which is the final output that contains 90 AAL Nodal efficiencies (one in each column) for all the subjects (one in each row).


