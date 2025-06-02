function e1 = GBS20P( )
%Batch testing script for xqsim program, tests +P method
%Matrix size       = 20*20
%Sequence length   = 1
%Total tests       = 9
%Total chi-2 error = 1.02
%Typical timing    = 5s

p.matrix     = @Identity;                        %matrix type
p.phase      = 1;
p.modes      = 20;                               %matrix size m
p.name       = sprintf('+P non-uniform GBS, M=%d',p.modes); 
p.tr         = 0.5*ones(1,p.modes);              %transmission
I            = ones(1,p.modes/5);                %identity vector
p.sqz        = [I/4,2*I,4*I,I/2,2*I];            %nonuniform squeezing
p.correl{4}  = 1:5;                              %correlation order 
p.correl{6}  = 1:5;                              %click order 
p.part{7}    = p.modes;                          %One-fold partititon 
p.part{8}    = p.modes/2*[1,1];                  %Two-fold partititon    
p.part{9}    = p.modes/4*[1,1,1,1];              %Four-fold partition
p.ensembles  = [1000,100,12];                    %repeats for errors
p.cutoff     = 1.e-7;
p.observe    = {@Np,@X2,@Y2,@Nm,@K,@Km,@K1,@Kn,@Kn};
p.compare    = {@Npc,@X2c,@Y2c,@Nmc,@Kc,@Kmc,@K1c,@Knc,@Knc};
p.glabels    = {{{},'Mode j'},{{},'Mode j'},{{},'Mode j'},{{},'Order'},...
               {{},'Mode j'},{{},'Clicks m'},{{},'Clicks m_1','Clicks m_2'},...
               {{},'Clicks m_1','Clicks m_2','Clicks m_3','Clicks m_4'}};
p.olabels    = {'<n>','<x^2>','<y^2>','<n_1...n_j>','<G^{(1)}>',...
               '<G^{(n)}>','G_1(m)','G_2(m)','G_4(m)'};
p.xk{4}      = p.correl(4);
p.xk{6}      = p.correl(6);
p.xk{7}      = {0:p.modes};
p.xk{7}      = {0:p.modes};
p.xk{8}      = {0:p.modes/2,0:p.modes/2};
p.xk{9}      = {0:p.modes/4,0:p.modes/4,0:p.modes/4,0:p.modes/4};
p0 = p;
[e1,d,cp]   = xqsim({p0});
xgraph(d,cp);
end