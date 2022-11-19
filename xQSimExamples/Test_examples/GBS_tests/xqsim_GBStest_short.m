function e1 = xqsim_GBStest_short( )
%Batch testing script for xqsim program
%Matrix size      = 20*20
%Sequence length  = 12
%Total tests      = 68
%Total graphs     = 152
%Chi-square error = 1.0 +/- 0.
%Typical timing   = 140s

p.matrix     = @Identity;                        %matrix type
p.M          = 20;                               %matrix size m
p.N          = 20;                               %input size n
p.name       = sprintf('+P non-uniform squeezed, M=%d',p.M); 
p.t          = 0.5;                              %transmission
I            = ones(1,p.M/5);                    %identity vector
p.r          = [I/4,2*I,4*I,I/2,I];              %nonuniform squeezing
p.eps        = zeros(1,p.M);                     %decoherence factor 
p.O{4}       = 1:5;                              %correlation order 
p.O{6}       = 1:5;                              %click order 
p.Part{7}    = p.M;                              %One-fold partititon 
p.Part{8}    = [p.M/2, p.M/2];                   %Two-fold partititon    
p.Part{9}    = [p.M/4, p.M/4,p.M/4, p.M/4];      %Four-fold partition
p.logs{9}    = [0,0,0,0,1];                      %Logs of data 
p.ensembles  = [100,10,12];                      %repeats for errors
p.cutoff     = 1.e-7;
p.observe    = {@n,@x2,@p2,@nm,@k,@km,@k1,@kn,@kn};
p.compare    = {@nc,@x2c,@p2c,@nmc,@kc,@kmc,@k1c,@knc,@knc};
p.glabels    = {{'Mode j'},{'Mode j'},{'Mode j'},{'Order'},...
               {'Mode j'},{'Clicks m'},{'Clicks m_1','Clicks m_2'},...
               {'Clicks m_1','Clicks m_2','Clicks m_3','Clicks m_4'}};
p.olabels    = {'<n>','<x^2>','<p^2>','<n_1...n_j>','<G^{(1)}>',...
               '<G^{(n)}>','G_1(m)','G_2(m)','G_4(m)'};
p.diffplot   = {2,2,2,2,2,2,2,2,2};
p.xk{4}      = p.O(4);
p.xk{6}      = p.O(6);
p.xk{7}      = {0:p.M};
p.xk{8}      = {0:p.M/2,0:p.M/2};
p.xk{9}      = {0:p.M/4,0:p.M/4,0:p.M/4,0:p.M/4};

p1           = p;
p1.eps       = 0.5*ones(1,p.M);                  %decoherence factor
p1.name      = sprintf('+P non-uniform thermalized, M=%d',p.M); %test name

p2           = p;
p2.r         = ones(1,p.M);
p2.eps       = zeros(1,p.M);                     %decoherence factor
p2.name      = sprintf('+P uniform squeezed, M=%d',p.M);     %test name

p3           = p;
p3.r         = ones(1,p.M);
p3.eps       = ones(1,p.M);                      %decoherence factor
p3.matrix    = @Unitary;                         %random unitary
p3.name      = sprintf('+P  unitary thermal, M=%d',p.M);      %test name

pW           = p;
pW.graphs    = 4;
pW.method    = 2; 
pW.name      = sprintf('W non-uniform squeezed, M=%d',p.M);  %test name

p1W          = p1;
p1W.graphs   = 4;
p1W.method   = 2;                              
p1W.name     = sprintf('W non-uniform thermalized, M=%d',p.M);%test name

p2W          = p2;
p2W.graphs   = 4;
p2W.method   = 2;                                 
p2W.name     = sprintf('W uniform squeezed, M=%d',p.M);      %test name

p3W          = p3;
p3W.graphs   = 4;
p3W.method   = 2;                              
p3W.name     = sprintf('W unitary thermal, M=%d',p.M);       %test name

pQ           = p;
pQ.graphs    = 4;
pQ.method    = 3; 
pQ.name      = sprintf('Q non-uniform squeezed, M=%d',p.M);  %test name

p1Q          = p1;
p1Q.graphs   = 4;
p1Q.method   = 3;                              
p1Q.name     = sprintf('Q non-uniform thermalized, M=%d',p.M);%test name

p2Q          = p2;
p2Q.graphs   = 4;
p2Q.method   = 3;                                 
p2Q.name     = sprintf('Q uniform squeezed, M=%d',p.M);      %test name

p3Q          = p3;
p3Q.graphs   = 4;
p3Q.method   = 3;                              
p3Q.name     = sprintf('Q unitary thermal, M=%d',p.M);       %test name

[e1,d,cp]    = xqsim({p,p1,p2,p3,pW,p1W,p2W,p3W,pQ,p1Q,p2Q,p3Q});
save xqsimtest.mat d cp
xgraph(d,cp);