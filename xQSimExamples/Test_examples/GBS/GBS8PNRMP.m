function e1 = GBS8PNRMP( )
%Batch testing script for xqsim PNR program using matrix-P
%Matrix size      = 8*8
%Sequence length  = 2
%Total tests      = 6
%Chi-square error = 1.27
%Typical timing   = 14s

p.matrix     = @Identity;                         %matrix type
p.modes      = 8;                                 %matrix size m
p.sqz        = ones(1,p.modes);                   %uniform squeezing
p.tr         = 1;
p.name       = sprintf('MP GBS PNR Ident M=%d tr=%d',p.modes,p.tr);
p.part{1}    = p.modes;                           %One-fold partititon 
p.part{2}    = p.modes/2*[1,1];                   %Two-fold partititon    
p.part{3}    = p.modes/4*[1,1,1,1];               %Four-fold partition
p.ensembles  = [10000,10,12];                     %repeats for errors
p.cutoff     = 1.e-5;
p.observe    = {@Mnn,@Mnn,@Mnn};                    
p.compare    = {@Nnsc,@Nnsc,@Nnsc};
p.glabels    = {{'Counts m'},{'Counts m_1','Counts m_2'},...
               {'Counts m_1','Counts m_2','Counts m_3','Counts m_4'}};
p.olabels    = {'G_1(m)','G_2(m)','G_4(m)'};
p.xk{1}      = {0:2*p.modes};
p.xk{2}      = {0:p.modes,0:p.modes};
p.xk{3}      = {0:p.modes/2,0:p.modes/2,0:p.modes/2,0:p.modes/2};
p1 = p;
p.tr         = 0.9;
p.compare    = {@Nnslc,@Nnslc,@Nnslc};
p.name       = sprintf('MP GBS PNR Ident M=%d tr=%d',p.modes,p.tr);
[e1,d,cp]    = xqsim(p1,p);
xgraph(d,cp);
end