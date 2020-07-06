function Zero_Diag=zero_diag(CorrMat)
SizeCorr=size(CorrMat,1);
CorrMat(1:SizeCorr+1:end)=0;
Zero_Diag=CorrMat;