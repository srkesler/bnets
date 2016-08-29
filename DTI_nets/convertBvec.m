% Shelli Kesler 8/20/16
% Convert FSL bvec files to Diffusion Toolkit format

[dirs] = textread('/Desktop/tracks/Paths.txt' , '%s');
% Txt file contains full paths to each subj e.g. /Volumes/STUDY/SUB001

for i = 1:length(dirs)
    cd(dirs{i});
    %5 = no. chars from end of the subj path that represent subj ID
    %starting count from 0, i.e. this selects the last 6 chars:
    bvec = importdata([dirs{i}(end-5:end) '_bvecs'],' ');  
    fileID = fopen('bvec.txt','w'); %create text file
    %format data in 3 columns separated by 1 space and a comma with 4 decimals
    fprintf(fileID,'%1.4f, %1.4f, %1.4f\n',bvec); 
end
