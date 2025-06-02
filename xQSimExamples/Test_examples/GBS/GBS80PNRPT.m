function e1 = GBS80PNRPT( )
%Batch testing script for xqsim program, +P, thermal input
%Matrix size      = 80*80
%Sequence length  = 1
%Total tests      = 6
%Chi-square error = 1.1 +/- 0.2
%Typical timing   = 23s

p.matrix     = @Identity;                        %matrix type
p.phase      = 1;
p.modes      = 80;                              %matrix size m
p.name       = sprintf('+P non-uniform GBS PNR M=%d',p.modes); 
p.tr         = 0.5*ones(1,p.modes);              %transmission
p.sqz        = ones(1,p.modes);                  %nonuniform squeezing
p.thermal    = ones(1,p.modes);                  %decoherence factor 
p.part{1}    = p.modes;                          %One-fold partititon 
p.part{2}    = p.modes/2;                        %One-fold partititon 
p.part{3}    = p.modes/2*[1,1];                  %Two-fold partititon 
p.part{4}    = p.modes/4*[1,1];                  %Two-fold partititon   
p.part{5}    = p.modes/4*[1,1,1];                %Three-fold partition
p.part{6}    = p.modes/5*[1,1,1,1];              %Four-fold partition
p.ensembles  = [1000,10,12];                     %repeats for errors
p.cutoff     = 1.e-6;
p.observe    = {@N1,@N1,@Nn,@Nn,@Nn,@Nn};
p.compare    = {@N1c,@N1c,@Nnc,@Nnc,@Nnc,@Nnc};
p.glabels    = {{{},'Counts m'},{{},'Counts m'},...
                {{},'Counts m_1','Counts m_2'},...
                {{},'Counts m_1','Counts m_2'}...
               {{},'Counts m_1','Counts m_2','Counts m_3'}...
               {{},'Counts m_1','Counts m_2','Counts m_3','Counts m_4'}};
p.olabels    = {'G_1(m)','G_1(m)','G_2(m)','G_2(m)','G_3(m)','G_4(m)'};
p.xk{1}      = {0:p.modes};
p.xk{2}      = {0:p.modes};
p.xk{3}      = {0:p.modes/2,0:p.modes/2};
p.xk{4}      = {0:p.modes/2,0:p.modes/2};
p.xk{5}      = {0:p.modes/4,0:p.modes/4,0:p.modes/4};
p.xk{6}      = {0:p.modes/4,0:p.modes/4,0:p.modes/4,0:p.modes/4};
p.diffplot   = {1,1,1,1,1,1};
p.name       = sprintf('matrix-P uniform GBS PNR M=%d',p.modes); 
[e1,d,cp]    = xqsim(p);
xgraph(d,cp);
end