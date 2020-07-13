function [pp1,c1,pp2,c2] = fdr_cor(p_value,ROI,Name,Pval_thr)
%input: accepts the vector of p-values (output of the regional graph analysis)(that comprises pos (G2>G1) and
%neg (G1>G2) p-vals) and computes the fdr corrected p-vals. It also acepts the
%matrix ROI which is a n*1 cell array containing the name of the ROIs. You
%need to choose a name for the measure you are analyzing (e.g
%'RegDegree_Group1')

%output: fdr_val: fdr corrected  p values

%%
% Hadi Hosseini (Jan 17, 2012)

%Pval_thr = 0.05; %for FDR (not for mafdr)

p1=p_value;p2=p_value;
for i=1:size(p_value,1)
    if p1(i) > 0
        p1(i) = p1(i) - 1;
    end
    if p2(i) < 0
        p2(i) = p2(i) + 1;
    end
end

%mafdar (does not work we1l) 
% [fdr1, fdr_p1] = mafdr(abs(p1), 'BHFDR', 'false', 'Lambda', [0.01:0.01:0.95], 'Method', 'bootstrap');
% [fdr2, fdr_p2] = mafdr(p2, 'BHFDR', 'false', 'Lambda', [0.01:0.01:0.95], 'Method', 'bootstrap');
% II = find(fdr_p1 < 0.05);
% JJ = find(fdr_p2 < 0.05);


% [c1,cc1] = FDR(abs(p1),Pval_thr);%c1(pos indep)  and cc1 (no assumption) are new thresholds based on FDR correction
% [c2,cc2] = FDR(p2,Pval_thr);

[h,c1,adj_p1] = fdr_bh(abs(p1),Pval_thr);
[h,c2,adj_p2] = fdr_bh(abs(p2),Pval_thr);


if ~isempty(c1)%(cc1)
    II = find(abs(p1) <= c1);
    % find regions with fdr vals < 0.05
    if ~isempty(II)
        fprintf('%-4s\n','');
        fprintf('%-4s\n',Name);
        fprintf('%-4s\n','');
        fprintf('%-4s\n',['group1 > group2; FDR corrected 0.05: ']); 
        fprintf('%-4s\n',''); 
        for i=1:size(II,1)
            fprintf('%-4s\n',cell2mat(ROI(II(i)))); 
        end
    end
end

if ~isempty(c2)%(cc2)
    JJ = find(p2 <= c2);
    if ~isempty(JJ)
        fprintf('%-4s\n','');
        fprintf('%-4s\n',Name);
        fprintf('%-4s\n','');
        fprintf('%-4s\n',['group2 > group1; FDR corrected 0.05: ']); 
        fprintf('%-4s\n','');
        for i=1:size(JJ,1)
            fprintf('%-4s\n',cell2mat(ROI(JJ(i)))); 
        end
        fprintf('%-4s\n','');
    end
end

%pp1 = abs(p1); pp2 = p2;
pp1 = adj_p1; pp2 = adj_p2;
save(['fdr_pval_Group1gdanGroup2_' Name], 'pp1', 'c1')
save(['fdr_pval_Group2gdanGroup1_' Name], 'pp2', 'c2')




