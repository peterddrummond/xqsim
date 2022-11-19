function e1 = xqsim_testcyc( )
%Cyclic testing script for xqsim program
%Matrix size      = 20*20
%Cycle length  = 40
%Total tests      = 1
%Total graphs     = 2
%Typical timing   = 10s

p.matrix     = @Identity;                        %matrix type
p.M          = 20;                               %matrix size m
p.N          = 20;                               %input size n
p.name       = sprintf('+P non-uniform squeezed, M=%d',p.M); 
p.t          = 0.7;                              %transmission
p.re         = 0.9;                              %recycling
p.cyc        = 20;                               %matrix size m
I            = ones(1,p.M/5);                    %identity vector
p.r          = [I/4,2*I,I,I/2,I];                %nonuniform squeezing
p.ensembles  = [1000,100,1];                      %repeats for errors
p.observe    = {@n};
p.glabels    = {{'cycle n','Mode j'}};
p.olabels    = {'<n>'};
[e1,d,cp]    = xqsim(p);
xgraph(d,cp);