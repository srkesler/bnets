function c_rand=hqs(c)

%Hirschberger-Qi-Steuer (HQS) Algorithm
%INPUT:
%c is the observed covariance matrix
%c=cov(randn(256,90));
%OUTPUT
%c_rand is a null correlation matrix for the observed covariance matrix. 

N=size(c,1); 

ind_upper=find(triu(ones(N,N),1)); 

e=mean(c(ind_upper));
if e<=0
    error('Off-diagonal mean less than or equal to zero');
end
v=var(c(ind_upper),1);
e_bar=mean(diag(c)); 
v_bar=var(diag(c),1);

m=max(2,round((e_bar^2-e^2)/v));
e_hat=sqrt(e/m); 
v_hat=-e_hat^2+sqrt(e_hat^4+v/m); 

Q=sqrt(v_hat)*randn(N,m)+e_hat;

c_rand=Q*Q'; 

%Convert to correlation matrices
c=diag(1./sqrt(diag(c)))*c*diag(1./sqrt(diag(c)));
c_rand=diag(1./sqrt(diag(c_rand)))*c_rand*diag(1./sqrt(diag(c_rand)));
