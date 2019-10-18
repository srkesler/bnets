% This script creates templates to extract cubes from a grey matter segmentation with dimension 91 x 109 x 91 voxels
%	The cubes have size 3 x 3 x 3 voxels, and therefore the template can start at different off set points.
%	This script provides separate templates for all off sets possible.
%
% Dependencies:
%	- standard matlab functions
%
% Adjust: Go to a directory where you have write access, where the bl_ind directory can live.
%
% Output : 	bl_ind directory
%		all_dirs.txt : overview of local directories --> corresponding to the possible off sets for that scans size
%	- bl_ind.m files for each off set --> array where each step of 27 indices = 1 cube.
%
% Author: Betty Tijms, 2011
% ----------------------------------------------------------------------------------------------------------------------------------

% Make a bl_ind directory
mkdir bl_ind_2mmiso
cd bl_ind_2mmiso

% Define the dimensions of the scans to be converted into networks
dimx=91;		
dimy=109;
dimz=91;

%define block size
n=3;

% Store all the directories in a text file --> needed when extracting network later
delete('all_dirs.txt');
fid1=fopen('all_dirs.txt', 'a');

% Loop through x, y and z dimensions to make sure to cover all possible combinations of [x, y, z] starting coordinates
for indx=1:n

	for indy=1:n

		for indz=1:n
			% Make the corresponding directory
			mkdir(strcat( char([int2str(indx), int2str(indy), int2str(indz)] ), '/data/rotation' ));
			mkdir(strcat( char([int2str(indx), int2str(indy), int2str(indz)] ), '/images' ));
			
			% Go to this directory
			cd( char([int2str(indx), int2str(indy), int2str(indz)]));

			store_dir= char([int2str(indx), int2str(indy), int2str(indz)]);
			fprintf(fid1, '%s \n', store_dir);

			%initialize variables
			b=zeros(dimx,dimy,dimz);	% Simulated volume with the same dimensions as the original volume
			tbl_ind=[];			% temporary cube index --> contains just one cube's worth of indices 
			bl_ind=0;			% the cube index, a concatenation of all the tbl_ind
			
			i_start=((indz-1)*(dimx*dimy))+((indy-1)*dimx)+indx;	%(note: indz -1 has been done, because the first z dimension is the first slice --> so no transitation is needed. Simarlaly for (indy-1)
			i_step=(dimx*dimy*n);
			end_i=numel(b)-i_step;
		
			% Determine end row 	
			row_n=floor(dimx/n);
			
			if (dimx-n*row_n+1)<indx
				end_k=n*row_n-n;	% --> for maximum 27 --> x max is 81
			else
				end_k=n*row_n;	% --> for maximum 27 --> x max is 81
			end
			%end_k=n*row_n-n;	% --> for 256-124 dimensions
		
			% Step size for columns: The step should encompase all the columns that belong to a cube
			j_step=n*dimx;
		
			% Determine end column --> this should be adjusted to the index where the cube gridding is starting.
			if (dimy-n*(floor(dimy/n) ) )< indy
				end_j=dimx*dimy-2*n*dimx;
			else
				end_j=dimx*dimy-n*dimx;
			end

			%data storage
			delete('data/bl_ind.m');
			fid2=fopen('data/bl_ind.m', 'a');
			
			% Loop through all the  voxels in all directions to extract cubes
			for i=i_start:i_step:end_i   
			
				for j=i:j_step:(i+end_j)
					for k=j:n:(j+end_k-1)
					% For each position on the z-axis extract a square of dimensions n x n
					
						for l=0:(n-1)
							% 0 is the first position on the z-axis
							tempz=dimx*dimy*l;
							
							start_k=k+tempz;
							% At this z-position, extract the square, by adding n-1 columns to this --> so 2:(n-1) times dimxz
							for m=0:(n-1)
								tempy=dimx*m;
								tbl_ind=[tbl_ind (start_k+tempy):(start_k+tempy+(n-1))];
							end
						end
					
					% store the indices corrsponding to cubes for this template
					fprintf(fid2, '%d \n', tbl_ind);
					%reset table_ind
					tbl_ind=[];
					%bl_ind=reshape(bl_ind',n,n,n);
					%b(tbl_ind)=k;
					end
				end
			end
			% close the file connection
			fclose(fid2);

			% Go back one directory
			cd ..
		end
	end
	indx
end
fclose(fid1)


% Batch script: Extracts networks from grey matter segmentations.
%	Network analyses is done with other scripts.
%
% Depends on:
% - the matlab nifti toolbox, download from: www/mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% - the output created with the 'create_cube_template.m script'.
% - SPM5 or SPM8 imcalc function
% - the functions that live in the "ALl_network_scripts" dir:
%	- determine_rois_with_minimum_nz : Determines which template extracts networks of minimum size.
%	- create_rrois : Permutes all the values in the scan.
%	- fast_cross_correlation: Computes the correlation matrices 
%
% PLEASE ADJUST the following lines:
% - line 43: add the path to the directory where you want to store the results
% - line 46: add the full path to the directory were All_network_scripts/functions/ live.
% - line 49: add the full path + file name of the txt file that contains the dir+names of all the scans (1 line per scan)
% - line 52: add the full path of the directory were the resliced images need to be written to.
% - line 55: add the full path of the directory where the spm canonical brain lives.
% - line 58: add the full path of the directory where /bl_ind directory lives, this directory contains the cube templates (output from 'create_cube_template.m')
%
%
%
% Output (stored in the subject directory unless specified otherwise):
% - iso2mm_s (subject number) .img : resliced images (this is stored in designated directory)
% - creates a directory for each subject (1:n). 
% - Sa = 3D volume that contains the MRI data in a matrix form
% - Va = header file that contains the scan info
% - nz = number of cubes (i.e., size of the network)
% - nancount = number of cubes that have standard deviation of 0. (these are excluded from further calculations)
% - off_set = the template used to extract the minimum number of cubes.
% - bind.mat (in single precision): this file contains the indices of the cubes that contain grey matter values
% - rois.mat and rrois.mat (in single precision): these files contain the cubes
% - lookup.mat and rlookup.mat : provide the link between the cubes from rois and rrois and the original scans (bind). These files are needed later when results are written in images (e.g., the degree or clustering values of the cubes).
% - rotcorr.mat and rrcorr.mat (in single precision) are the correlation matrices maximised for with reflection over all angles and rotation with multiples of 45degrees.
% - th = threshold to binarise the matrices: corresponds to the correlation value where in the random image 
% - fp = percentage of positive random correlations
% - sp = sparsity of the matrix with threshold p_corrected = 0.05
%
% Author: Betty Tijms, 2011, version 12.08.2013
% -------------------------------------------------------------------------------------------------------------------------------

%cd /usr/local/MATLAB/R2011a/bin/
% matlab -nodesktop

% --- ADJUST FOLLOWING ----%
result_dir ='/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/results_2mmiso';
mkdir(result_dir)
cd(result_dir)
spm8
% FILL IN THE FILE NAME with full path that contains the text file with all the grey matter segmentations 
[CN_a1]=textread('/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/filelist_GM.txt','%s');


% FILL IN THE FULL PATH of the directory where the resliced images should live.
reslice_dir = '/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/reslice_2mmiso/';
mkdir(reslice_dir)
% Add the path to the Nifti toolbox (see read me file for download link)
addpath /Volumes/wdext/scripts/matlab_shared_toolboxes/NIfTI_20140122/

% PATH where the canonical SPM image "avg305T1.nii´ can be found --> this is used to reslice the images
P1= '/Volumes/wdext/scripts/matlab_shared_toolboxes/spm8/apriori/brainmask.nii';

% Please provide path to the folder Extract_individual_GM_networks
all_network_scripts_path ='/Volumes/wdext/scripts/Single_Subject_Grey_Matter_Networks-master/Extract_individual_GM_networks_v20150902/';

% Please provide path to the folder bl_ind
bl_dir = '/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/bl_ind_2mmiso/';


%% -----Finished with the adjustments ---%

% add path to the functions in All_network_scripts folder
addpath(strcat(all_network_scripts_path,'functions'))

% Get content of all_dirs.txt files ( this file contains all the dirs of the cube template for all possible off set values)
CN_a2 =textread(strcat(bl_dir, 'all_dirs.txt'),'%s');

% end of loop
numimages=size(CN_a1,1);

% number of dimensions
n=3;
s=n^3;

% Loop through all the grey matter segmentations listed in CN_a1 and extract network. Takes ~25 minutes per scan.

for im=1:numimages	
	% Get this scan and do the following loop to threshold and then reslice it
	this_scan=char(CN_a1(im));
	
	% make the directory for this subject
	mkdir(strcat('s', int2str(im),'/data/rotation/'));
	mkdir(strcat('s', int2str(im),'/images/'));
	
	% Go to this directory
	cd(strcat('s', int2str(im)));
	
	% make a temporary copy of the grey matter segmentation in this directory
	copyfile(this_scan, 'temp_grey.nii');
	
	% which directory are we now?
	t_result_dir = strcat(result_dir, '/s',int2str(im));
% % % 	
% % % 	% Make sure that the native image realigns with MNI images
% % % 	coreg_est = spm_coreg(P1, t1_scan);
% % % 	M = inv(spm_matrix(coreg_est));
% % % 	MM = spm_get_space('temp_grey.nii'); %adjust the header of the grey matter image to store these realigned parameters
% % % 	spm_get_space('temp_grey.nii', M *MM);
	
% % % 	% save the reorientation paramters
% % % 	save images/coreg_est.mat coreg_est
	
%  	%% Next reslice the scan to 2x2x2 mm isotropic voxels to reduce amount of data.
  	P = strvcat(P1, 'temp_grey.nii');
  	Q = strcat(reslice_dir,'/iso2mm_s', num2str(im), '.nii');
  	f = 'i2';
  	flags = {[],[],[],[]};
  	R = spm_imcalc_ui(P,Q,f,flags);
%  
%  	%close the figure window
 	close
 	
 	% remove the temporary grey matter image
 	delete('temp_grey.nii');

%	%% Now extract Sa and Va --> Va contains all the info from the .hrd file, Sa is the actual image, it contains the grey matter intensity values.
	Va=spm_vol(strcat(reslice_dir, '/iso2mm_s', num2str(im), '.nii'));
% % %     Va=spm_vol(this_scan); % added this in the modified version by Vikram on 12/07/15
	Sa=spm_read_vols(Va);

	% Save them in the current directory
	save Va.mat Va;
	save Sa.mat Sa;

	% clear variables that aren't needed anymore
	clear this_scan P Q R f flags;

	%% Get the off_set bl_ind (this corresponds to the indices that make up the 3d cubes)  with the minimum number of cubes (min. nz).
	tic
    [nz, nan_count, off_set] = determine_rois_with_minimumNZ(CN_a2, bl_dir, n, Sa, t_result_dir);
toc
	% Store nz.m : the number of cubes, this is the size of the network
	save data/nz.mat nz;

	% nan_count = number of cubes that have a variance of 0, these are excluded because cannot compute correlation coefficient for these cubes.
	save data/nan_count.mat nan_count;
	% off_set indicates which specific template was used to get the cubes.
	save data/off_set.mat off_set;

	% Get bind and store bind too. Bind contains the indices of the cubes, so we can efficiently do computations later. It is a long vector of which every 27 consecutive voxels are a cube.
	%tt=strcat(bl_dir, off_set, '/data/bind.m');
	load data/bind.m;
	
	% Convert to single, for memory reasons
	bind=single(bind);
	delete data/bind.m;
	save data/bind.mat bind	;

	% Now create:
	% - the rois: This is a matrix of which each column corresponds to a cube, we will compute use this to compute the correlations
	% - lookup table --> to lookup the cubes that correspond to the correlations.
	
	lb=length(bind); % End point for loop
	col=1;		% iteration counter
		
	rois=zeros(s,nz, 'single');	% Create variable to store the rois (i.e., cubes)
	lookup=zeros(s,nz, 'single');	% Create variable to store the lookup table --> this table links the cubes to the indices in the original iso2mm image. 
	
	% The next loop gos through bind for indices of voxels that belong to each ROI (i.e. cube)
	% lb is the last voxel that belongs to a roi
	%lookup is lookup table to go from corr index to bind to Sa
	
	for i=1:s:lb	
		rois(:,col)=Sa(bind(i:(i+(s-1))));
		lookup(:,col)=i:(i+(s-1));
		col=col+1;
	end
	

	% Save the rois and lookup table
	save data/rois.mat rois;
	% look up table --> to go from corr index to bind to Sa
	save data/lookup.mat lookup;

	%clear unneeded variables
	clear lookup lb ;

	%% Randomise ROIS: Create a 'random brain' to estimate the threshold with for later stages
	[rrois, rlookup]= create_rrois (rois, n, Sa, Va, off_set, bind, bl_dir, nz);

	% save rrois and rlookup
	save data/rrois.mat rrois;
	save data/rlookup.mat rlookup;

	% Remove all variables that take space and aren't needed anymore
	clear bind lookup rlookup Sa Va;
	
	% Set the following variable to 2 when correlation is maximised for rotation with angle multiples of 45degrees.
	forty=2;
%     parpool;
    tic
	[rotcorr, rrcorr] = fast_cross_correlation(rois, rrois, n,forty);
	toc
	% Save it
	save data/rotation/rotcorr.mat rotcorr;
	save data/rotation/rrcorr.mat rrcorr;

	%% Now get the threshold and save it
	[th, fp, sp] = auto_threshold(rotcorr, rrcorr, nz);
    rotcorr_th=single(rotcorr>th);
    % Add this threshold to all_th
	%all_th(im,:)=[th,fp,sp];
	%save this and get later
	save data/th.mat th;
	save data/fp.mat fp;
	save data/sp.mat sp;
    save data/rotation/rotcorr_th.mat rotcorr_th;

	% Clear all variables
	clear rotcorr rotcorr_th rrcorr th fp sp ;


	cd ..
	im

end




% Link cubes to AAL mask --> to be able to extract % hubs e.g.


% set result dir
result_dir = '/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/results_2mmiso';
% reslice_dir = '/Volumes/wdext/Study/Glioma/New/single_subject_gm_network/results_2mmiso/';
roi_folder='/Volumes/wdext/scripts/rois/aalroi_for_tijmstoolbox/91x109x91_orig/'; % no need to change this. 
% Please enter the number of subjects to process below
numimages= 71; %105

% im= 134 = excluded
for im=1:numimages
	% Go to this person directory
	cd(strcat(result_dir, '/s',int2str(im)));
	%mkdir atlas

	% Get the mask : nifti toolbox does not always work properly -> in that case use SPM load images.
	%taSa = load_nii(strcat(reslice_dir, 'iso2mm_atlas_', char(all_ids(im)),'.nii' )); 
	load(strcat(roi_folder,'aal_90.mat'));	% contains the image matrix of the resliced atlas image, as save after loading with SPM. 
	%aSa=aSa.img;
	aal = round(aal);	% make sure that the images contain only round numbers from 1:90

	% Load other data necessary from this subject
	load Va.mat
	load data/bind.mat
	load data/lookup.mat
	load data/nz.mat
	
	% Label the voxels in a roi
	label_rois = zeros(size(lookup));

	s = 27;
	lb=length(bind); % End point for loop
	col=1;		% iteration counter
	
	for i=1:s:lb	
		label_rois(:,col)=aal(bind(i:(i+(s-1))));
		col=col+1;
	end
			
	% Remove cubes that are not labelled 
	aal_rois_ind = find(sum(label_rois)~=0); % removes zero-labelled cubes
	aal_rois_labels = label_rois(:,aal_rois_ind);
	aal_rois_labels (aal_rois_labels ==0)=NaN;
	aal_rois_labels = mode(aal_rois_labels); %finds the label that has max number (mode) of voxels within a cube
	
	% Save indices of labeled rois --> needed because we only want the stats from this subset.	
	save aal_rois_ind.mat aal_rois_ind   % This is the index of nonzero label cubes. For exampleout of 8668 only 6939 cubes might have label. This is the index of which ones are nonzero cubes.
	save aal_rois_labels.mat aal_rois_labels % this is the label of nonzero cubes. each column will have a number in the range 1-90. 
    % The label names can be accessed from: /Volumes/wdext/scripts/rois/aalroi_for_tijmstoolbox/91x109x91_orig/aal_label_name.mat
	cd ..

	im
end

% Shelli Kesler and Vikram Rao; created on  8/22/16
% Conducts targeted network attack by iteratively removing nodes in
% descending order of degree. Global efficiency is calculated after each
% attack iteration and then area under the curve is computed across all
% iterations. Here we have used AAL regions from Tijms' toolbox. We start
% with say about 8000 nodes and then each cube (node) is assigned an AAL
% area. There might be some nodes for which there are no AAL regions
% assigned, such nodes are first removed. Next on the remaining nodes
% degrees are calculated for each node. Next degrees of similar AAL cubes
% are averaged. For example, for all nodes that are assigned to Left
% Hippocampus, an mean degree is calculated. We end up with an array of
% 90 that has mean degree for each roi. Next they are arranged in
% descending order of degree. In each iteration (total 90) henceforth, all nodes that
% belong to the highest mean degree roi are removed (zeroed) and global
% efficiency is calculated, next the second highest mean degree along
% with the highest mean degree roi nodes are removed (zeroed) and global
% efficiency is again calculated for that iteration. After all iterations
% are done, a AUC is calculated for that subject. Variables saved:
% 1. TAgeff: Global efficiency after each attack. This will be a 90 X 1 array
% 2. roi_ind: This is the index of each roi (1 to 90) in descending order of mean degree 
% 3. aal_roimean_deg_sorted: This is the actual mean degree for each roi sorted descendingly
% 4. TAgeffAUC_subj: This is the AUC for each subject 
% 5. rotcorr_th_labeled: This is the version of rotcorr_th matrix after removing the unlabeled nodes

clear all;
clc;
% Location to store the final AUC for all subjects
root='/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/results_2mmiso';
% Adding paths
addpath(genpath('/Volumes/wdext/scripts/matlab_shared_toolboxes/BCT/2016_01_16_BCT/'))
spm8
% Location to text file where the paths for each subject's Single subject
% GM network analysis data are stored
[subjs]=textread('/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/results_2mmiso/paths.txt','%s');
%% This is the main loop that loops across each subject
for j = 1:length(subjs)
    j
    cd (fullfile(subjs{j,1},'data')); % Cd into each subjects directory
    disp(sprintf('working on subject s%d',j))
    % load the correlation matrices along with the aal assignment for each node
    load(fullfile(subjs{j,1},'data','rotation','rotcorr_th.mat'));
    load(fullfile(subjs{j,1},'aal_rois_ind.mat'));
    load(fullfile(subjs{j,1},'aal_rois_labels.mat'));
    %Remove the unlabeled nodes from the rotcorr_th matrix and call it 'rotcorr_th_labeled'
    rmvnonlblnodes=sort(setdiff(1:length(rotcorr_th),aal_rois_ind),'descend'); %Finds indices of unlabeled nodes
    rotcorr_th_labeled=rotcorr_th;
    rotcorr_th_labeled(rmvnonlblnodes,:)=[]; rotcorr_th_labeled(:,rmvnonlblnodes)=[]; %Removes those columns and rows starting from the end
    %Calculate degree of each labeled nodes 
    deg_labeled= degrees_und(rotcorr_th_labeled);
    sorted_aal_rois_labels=unique(aal_rois_labels); % Sorts labels from 1 to 90
    %Calculate mean degree of all nodes in each roi
    for nrois=1:length(unique(aal_rois_labels))
        iroi_ind=find(aal_rois_labels==sorted_aal_rois_labels(nrois)); % Finds nodes indices for each roi(label)
        aal_roimean_deg(nrois,1)=mean(deg_labeled(iroi_ind));
    end
    %Arrange the mean degree of each roi in descending order
%     [aal_roimean_deg_sorted,roi_ind]=sort(aal_roimean_deg,'descend');
    
    %Initialize the matrix 'rotcorr_th_labeled' by making a copy of it first
    G = rotcorr_th_labeled;
    tic; localeff = efficiency_bin(G,1); toc;
    
    for nrois=1:length(unique(aal_rois_labels))
        iroi_ind=find(aal_rois_labels==sorted_aal_rois_labels(nrois)); % Finds nodes indices for each roi(label)
        aal_roimean_localeff(nrois,1)=mean(localeff(iroi_ind));
    end
% %     for irois=1:length(roi_ind)
% %         iroi_ind=find(aal_rois_labels==sorted_aal_rois_labels(roi_ind(irois))); % Finds nodes indices for each roi(label)
% %         G(iroi_ind,:) = 0; G(:,iroi_ind) = 0; % zero out the nodes associated with each roi
% %         TAgeff(:,irois) = efficiency_bin(G,1); % Calculate global efficiency
% %     end
    %Loops across each roi (90) in the descending order of mean degree and
    %zeros out the nodes associated with each roi and then find the global
    %efficiency for that roi. Does this iteratively until the last roi.
    %Note: For example in the 5th iteration, all nodes associated with the top 5 rois (descending) in order of mean degree are zeroed out.  
    
    save AAL_localeff.mat localeff aal_roimean_deg aal_roimean_localeff;%save variables of interest
    clear  aal_roimean_deg aal_roimean_localeff G % Clear variables after each iteration   
end


% Shelli Kesler 2/2/2016

% Extracts in a single.mat file nodal efficiencies of all 90 rois for all
% the subjects in the study

clear all;
clc;
% Location to store the final AUC for all subjects
root='/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/results_2mmiso';
% Adding paths
addpath(genpath('/Volumes/wdext/scripts/matlab_shared_toolboxes/BCT/2016_01_16_BCT/'))
spm8
% Location to text file where the paths for each subject's Single subject
% GM network analysis data are stored
[subjs]=textread('/Volumes/wdext/Study/Glioma/Glioma71/single_subject_gm_network/results_2mmiso/paths.txt','%s');
for j = 1:length(subjs)
    load(fullfile(subjs{j,1},'data','AAL_localeff.mat'));
    nodal_efficiency(j,1:90)=aal_roimean_localeff';
    meandeg(j,1)=mean(aal_roimean_deg);
    nodes(j,1)=length(localeff);
end
nodalefficiency=[nodal_efficiency meandeg nodes]; % Combining degree, size and nodal efficiency
save(fullfile(root,'nodalefficiency_stats.mat'), 'nodalefficiency');       




