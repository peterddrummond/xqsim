function e1 = xqsim_HOC_test20( )
%Testing script for xqsim higher order correlations program
%Matrix size      = 20*20
%Sequence length  = 1
%Total tests      = 5
%Total graphs     = 10
%Chi-square error = 1.0 +/- 0.1
%Typical timing   = 76s

p.matrix     = @Unitary;                         %matrix type
p.M          = 20;                               %matrix size m
p.N          = 20;                               %input size n
p.name       = sprintf('+P  unitary thermal, M=%d',p.M); 
p.t          = 0.5;                              %transmission
I            = ones(1,p.M/5);                    %identity vector
p.r          = ones(1,p.M);                      %nonuniform squeezing
p.eps        = ones(1,p.M);                      %decoherence factor 
p.CO{2}      = 2;                              %correlation order 
p.CO{3}      = 3;                              %click order 
p.CO{4}      = 2;
p.Part{5}    = p.M;                              %One-fold partititon 
p.ensembles  = [1000,100,12];                    %repeats for errors
p.cutoff     = 1.e-7;
p.observe    = {@k,@km2,@km3,@kmsub,@k1};
p.compare    = {@kc,@km2c,@km3c,@kmsubc,@k1c};
p.glabels    = {{'j'},{'j<k'},{'j<k<h'},{'j,k'},{'Clicks m'}};
p.olabels    = {'<\pi_j(1)>','<\pi_j(1)\pi_k(1)>','<\pi_j(1)\pi_k(1)\pi_h(1)>',...
                '<\pi_j(1)\pi_k(1)>', 'G_1(m)'};
p.diffplot   = {2,2,2,2,2};
p.xk{2}      = {0:190};
p.xk{3}      = {0:1140};
p.xk{4}      = {0:p.M-1};
p.xk{5}      = {0:p.M};
[e1,d,cp]    = xqsim(p);
%save xqsimtest.mat d cp
xgraph(d,cp);