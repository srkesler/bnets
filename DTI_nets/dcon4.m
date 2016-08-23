clear all;

% Shelli Kesler 6/19/15
% Calculates mean volume for each pair of ROIs, imports 
% UMCP fiber number (FN) and stat txt (FA) files, thresholds the fiber 
% count, weights the fiber count by the stat values and corrects these 
% products for average ROI volume.

% This script assumes that dcon1, dcon2 and dcon3 have been 
% completed and that subject folders are located in group specific folders.

% This script requires the txt2mat script in the Matlab path
% http://www.mathworks.com/matlabcentral/fileexchange/18430-txt2mat

% Shelli Kesler 8/28/15
% Updated to use fslstats ROI volume output

% paths to subject folders
[subjs] = textread('fullpath/paths.txt','s%');

%path to folder that will contain all FNFA*.mats
datadir2 = '/path/FNFA/'; 
mkdir(datadir2);

% calculate mean volume for each pair of ROIs
for j=1:length(subjs)
    cd(subjs{j});
    c = csvread('ROIvols.csv'); %output from fslstats
    r = c';
    n = 1;
    i=1:90;
    while n <= 90
        ROIvols(i,n)=(c(i)+r(:,n))/2;
        n = n+1;
    end;
    save ROIvols.mat ROIvols;
    
    % convert subjconnectmat.txt (FN) to mat file and save it
    subjconnectmat=txt2mat('subj_connectmat.txt');
    save subjconnectmat.mat subjconnectmat;

    % set FN variable in Matlab environment
    subjconnectmat_thresh3 = subjconnectmat;
    
    % fiber number false positive threshold
    subjconnectmat_thresh3(subjconnectmat < 3) = 0; 
    
    % save thresholded FN matrix
    save subjconnectmat_thresh3.mat subjconnectmat_thresh3;

    % convert subjstatmat.txt (FA) to mat file and save it
    subjstatmat=txt2mat('subj_statmat.txt');
    save subjstatmat.mat subjstatmat;

    % multiply FN by FA and divide by ROI volume
    FNFAtemp = times(subjconnectmat_thresh3,subjstatmat);
    FNFA = FNFAtemp./ROIvols;
    
    % save all FNFA.mat's in one folder
    subjno = num2str(j,'%03d');
    save(['FNFA_' subjno], 'FNFA');
    copyfile(['FNFA_' subjno '.mat'], datadir2);
end;
