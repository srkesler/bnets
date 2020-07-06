function regMes_long_1g(data_log,Data1,Data2)

% Edit here------------------------

ff = load(data_log,'mat4GAT');
G1T1_temp=load(Data1,'NetMes_Bin'); G1T1=G1T1_temp.NetMes_Bin;
G1T2_temp=load(Data2,'NetMes_Bin'); G1T2=G1T2_temp.NetMes_Bin;
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
cd Regional
save RegMesAUC EdgeAUCG1T1 EdgeAUCG1T2