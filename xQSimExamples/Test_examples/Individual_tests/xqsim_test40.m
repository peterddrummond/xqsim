function e1 = xqsim_test40( )
%Testing script for xqsim program
%Matrix size      = 40*40
%Sequence length  = 1
%Total tests      = 8
%Total graphs     = 18
%Chi-square error = 1.0 +/- 0.1
%Typical timing   = 6s

p.matrix     = @Identity;                        %matrix type
p.M          = 40;                               %matrix size m
p.N          = 40;                               %input size n
p.name       = sprintf('+P non-uniform squeezed, M=%d',p.M); 
p.t          = 0.5;                              %transmission
I            = ones(1,p.M/5);                    %identity vector
p.r          = [I/4,2*I,4*I,I/2,I];              %nonuniform squeezing
p.eps        = zeros(1,p.M);                     %decoherence factor 
p.O{4}       = 1:5;                              %correlation order 
p.O{6}       = 1:5;                              %click order 
p.Part{7}    = p.M;                              %One-fold partititon 
p.Part{8}    = [p.M/2, p.M/2];                   %Two-fold partititon    
p.ensembles  = [1000,100,12];                    %repeats for errors
p.cutoff     = 1.e-7;
p.observe    = {@n,@x2,@p2,@nm,@k,@km,@k1,@kn};
p.compare    = {@nc,@x2c,@p2c,@nmc,@kc,@kmc,@k1c,@knc};
p.glabels    = {{'Mode j'},{'Mode j'},{'Mode j'},{'Order'},...
               {'Mode j'},{'Clicks m'},{'Clicks m_1','Clicks m_2'}};
p.olabels    = {'<n>','<x^2>','<p^2>','<n_1...n_j>','<G^{(1)}>',...
               '<G^{(n)}>','G_1(m)','G_2(m)'};
p.diffplot   = {2,2,2,2,2,2,2,2};
p.xk{4}      = p.O(4);
p.xk{6}      = p.O(6);
p.xk{7}      = {0:p.M};
p.xk{8}      = {0:p.M/2,0:p.M/2};

[e1,d,cp]    = xqsim(p);
save xqsimtest.mat d cp
xgraph(d,cp);