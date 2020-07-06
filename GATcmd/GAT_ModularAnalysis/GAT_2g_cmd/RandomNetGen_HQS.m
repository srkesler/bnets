function r_null_unthr = RandomNetGen_HQS (CX,Iter)

nROI = size(CX,1);



parfor i = 1:Iter

    r_null = hqs (CX);
    
    
    r_null(1:nROI+1:end) = 0;

    r_null(r_null<0) = 0;
    
    R_null(:,:,i) = r_null;

end
r_null_unthr = R_null;

