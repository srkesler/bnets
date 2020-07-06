function regMes_2g(data_log,Data1,Data2,outname)

% Edit here------------------------

ff = load(data_log,'mat4GAT');
G1_temp=load(Data1,'NetMes_Bin'); G1=G1_temp.NetMes_Bin;
G2_temp=load(Data2,'NetMes_Bin'); G2=G2_temp.NetMes_Bin;
%-----------------------------------

for sub = 1:size(G1,1)
    
    for i=1:size(G1,2) %density loop
        EdgeG1(:,:,i) = G1{sub,i}{18,3};
    end
    EdgeAUCG1(:,:,sub) = trapz(EdgeG1,3);
end

for sub = 1:size(G2,1)
    
    for i=1:size(G2,2) %density loop
        EdgeG2(:,:,i) = G2{sub,i}{18,3};
    end
    EdgeAUCG2(:,:,sub) = trapz(EdgeG2,3);
end
cd Regional
save(sprintf('RegMesAUC_%s',outname), 'EdgeAUCG1', 'EdgeAUCG2','-v7.3')