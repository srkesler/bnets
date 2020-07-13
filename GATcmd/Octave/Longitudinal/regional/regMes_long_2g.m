function regMes_long_2g(data_log,Data1_g1g2,Data2_g1g2,Data3_g3g4,Data4_g3g4)

% Edit here------------------------

ff = load(data_log,'mat4GAT');
G1T1_temp=load(Data1_g1g2,'NetMes_Bin'); G1T1=G1T1_temp.NetMes_Bin;
G1T2_temp=load(Data2_g1g2,'NetMes_Bin'); G1T2=G1T2_temp.NetMes_Bin;
G2T1_temp=load(Data3_g3g4,'NetMes_Bin'); G2T1=G2T1_temp.NetMes_Bin;
G2T2_temp=load(Data4_g3g4,'NetMes_Bin'); G2T2=G2T2_temp.NetMes_Bin;
%-----------------------------------

for sub = 1:size(G1T1,1)
    
    for i=1:size(G1T1,2) %density loop
        EdgeG1T1(:,:,i) = G1T1{sub,i}{18,3};
    end
    EdgeAUCG1T1(:,:,sub) = trapz(EdgeG1T1,3);
end

for sub = 1:size(G1T2,1)
    for i=1:size(G1T2,2) %density loop
        EdgeG1T2(:,:,i) = G1T2{sub,i}{18,3};
    end
    EdgeAUCG1T2(:,:,sub) = trapz(EdgeG1T2,3);
end

for sub = 1:size(G2T1,1)
    for i=1:size(G2T1,2) %density loop
        EdgeG2T1(:,:,i) = G2T1{sub,i}{18,3};
    end
    EdgeAUCG2T1(:,:,sub) = trapz(EdgeG2T1,3);
end

for sub = 1:size(G2T2,1)
    for i=1:size(G2T2,2) %density loop
        EdgeG2T2(:,:,i) = G2T2{sub,i}{18,3};
    end
    EdgeAUCG2T2(:,:,sub) = trapz(EdgeG2T2,3);
end

cd Regional
save RegMesAUC.mat -v7 EdgeAUCG1T1 EdgeAUCG1T2 EdgeAUCG2T1 EdgeAUCG2T2