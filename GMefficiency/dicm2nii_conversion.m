clc
clear all
close all
spm8
addpath('/Volumes/wdext/scripts/matlab_shared_toolboxes/dicm2nii');
filelist_name='/Volumes/wdext/Study/Glioma/Glioma71/scripts/filelist.txt';
fid = fopen(filelist_name);
file_list = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
file_list=file_list{1,1};
for idir=1:71
    idir
    dicm2nii(fullfile(file_list{idir},'anat'), fullfile(file_list{idir},'anat'),'nii');
    cd(fullfile(file_list{idir},'anat'))
    h=dir('*.nii');
    movefile(h.name,'anat.nii')
end
