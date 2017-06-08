function CL=CL_per(a)

%a:n*m with n permutations for m variables (nets) 

[n,m]=size(a);
lim=round(.05*n);
CL=[];
for i=1:m
    
    xx = a(:,i);
    xx(isnan(xx)) = nanmean(xx);
    MM = mean(xx);
    xx = xx - MM;
   
    
    [as_abs,ix] = sort(abs(xx));
    
    ix_ac = ix(n:-1:n-lim-1);
    ac = xx(ix_ac);
    
    limUp = size(find(ac >= 0),1);
    limLo = size(find(ac < 0),1);
    
    as = sort(xx);
    
    if n-limUp+1 <= n
        CLup=as(n-limUp+1);
    else
        CLup=as(n);
    end
    
    CLlo=as(limLo+1);
    
    CL=[CL;CLup+MM CLlo+MM];
    
end
CL=CL';
