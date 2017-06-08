% Shelli Kesler

clear all;
clc;

% Measure small-world properties across various densities from minimum to
% maximum connection density (or user specified values) and compare the AUCS between groups using
% permutation testing

addpath('/Users/skesler/Dropbox/MATLAB/myScripts/bNets/swiAUC'); %path to dependency functions
addpath('/Users/skesler/Dropbox/MATLAB/BCT/2017_01_15_BCT') % path to brain connectivity toolbox
tic
load('Z8groupsTgNtg.mat'); % EDIT: concatenated, square matrices NOTE: matrices must be ordered by group and labeled Z
densStep = .02; % EDIT: density interval step
Nperm = 1000; % EDIT: number of random bootstrap networks
G1N = 4; % EDIT: number of subjects in group 1
swiLabels = {'swi';'cc';'pl';'geff';'leff';'mod'}; % EDIT based on any additions to swiFx

% find minimum connection density across subjects
[minD,V1,D] = minDensity(Z);
cMinD = minD;
minD = 0.1; % EDIT: comment out or enter desired minimum density value

%Find maximum density
for i = 1:size(Z,3)
    X = Z(:,:,i);X(X<0)=0;X(logical(eye(size(X)))) = 0;
    iDens(i,:) = density_und(X);
end
maxD = max(iDens);
cMaxD = maxD;
maxD = 0.8; %EDIT: comment out or enter desired max density value

steps = minD:densStep:maxD;
NSUB = size(Z,3);
NSTEP = length(steps);
ZTmats = cell(NSUB,NSTEP);
% swiData = cell(NSUB,NSTEP);
swiData = cell(1,NSTEP); % changed by VR on 06/07/2017
isempty(gcp);

% calculate connectome properties across densities
parfor interval=1:NSTEP
  S = steps(interval);
  ZT = binGraph(Z,V1,D,S);
  ZTmats{interval} = ZT; %ZT = 1xnsub cell where ZT{1,subj} = subj's thresholded matrix for S density
  swi = swiFx(ZT);
  swiData{interval} = swi; %swiData{1,density}(:,measure) = measures for each density for all subjs; order of measures in swiLabels
end

% calculate group difference in AUCs
for a = 1:NSTEP
    for m = 1:size(swiData{1,1},2)
    swiMean1(m,a) = mean(swiData{1,a}(1:G1N,m));
    swiMean2(m,a) = mean(swiData{1,a}((G1N+1:end),m));
    swiAUCdiff(m,:) = abs(trapz(swiMean1(m,:))-trapz(swiMean2(m,:))); %nMeasurex1 vector
    end
end

% generate random boostrapped networks
for rnets = 1:Nperm
    randG{rnets} = datasample(Z,NSUB,3,'Replace',true); %1xNperm cell where randG{1,Nperm} contains nxnxm Z matrix for a random group
end

% get densities for random networks
parfor randperm = 1:Nperm
    Zrand = randG{1,randperm}; %nxnxm Z matrix from randperm permutation
    [grpDens,V1rand] = randDensity(Zrand,NSUB);
    DensR{randperm} = grpDens; %1xNperm cell where DensR{1,Nperm}= Nperm random group's densities
    V1R{randperm} = V1rand; %1xNperm cell where V1R{1,Nperm}=Nperm random group's sorted connectivity values
end

% measure connectome properties for bootstrapped networks
parfor nrand=1:Nperm
    Zrand2 = randG{1,nrand};
    for intR = 1:NSTEP
        Srand = steps(intR);
        ZTrand = binRandGraphs(Zrand2,V1R,DensR,Srand,Nperm);
        ZTmatsRand{nrand,intR} = ZTrand; %npermXsteps cell where {nperm,step} = 1xnsub cell with matrices for each subj at step
        swiRand = swiFx(ZTrand);
        swiDataRand{nrand,intR} = swiRand; %npermXsteps cell where {nperm,step} = nsub x nMeasure swiData for nperm random group
    end
end

% calculate random group differences in AUCs
for nperm = 1:Nperm
    for wStep = 1:NSTEP
        for mes = 1:size(swiDataRand{1,1},2)
            swiRmean1(mes,wStep) = mean(swiDataRand{nperm,wStep}(1:G1N,mes));
            swiRmean2(mes,wStep) = mean(swiDataRand{nperm,wStep}((G1N+1:end),mes));
            randAUCdiff(mes,:) = abs(trapz(swiRmean1(mes,:))-trapz(swiRmean2(mes,:)));
        end
    end
    swiAUCdiffRand(nperm,:) = randAUCdiff; %nperm x nMeasure matrix i.e. each column is permutation distribution for that measure
end

% calculate pvalues and confidence intervals
lolim = round(.025*Nperm); hilim = round(Nperm-lolim);
for iMes = 1:size(swiDataRand{1,1},2)
    permPvals(iMes,:) = mean(swiAUCdiffRand(:,iMes) >= swiAUCdiff(iMes,:)); %two-tailed
    Mtemp=sort(swiAUCdiffRand(:,iMes)); 
    CIs(1,iMes)=Mtemp(lolim); CIs(2,iMes)=Mtemp(hilim+1);
end
CI_gat = CL_per(swiAUCdiffRand);

save swiAUC.mat ZTmats swiData swiLabels swiMean1 swiMean2 swiAUCdiff minD...
    maxD randG DensR ZTmatsRand swiDataRand swiAUCdiffRand permPvals CIs CI_gat...
    cMinD cMaxD -v7.3;
toc