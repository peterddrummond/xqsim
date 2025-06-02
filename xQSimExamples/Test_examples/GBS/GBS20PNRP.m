function e1 = GBS20PNRP( )
%Batch testing script for xqsim program, PNR pos P tests
%GBS PNR example with 50% transmission
%Matrix size       = 20*20
%Sequence length   = 1
%Total tests       = 1
%Total chi-2/k     = 0.27
%Typical timing    = 7s

p.matrix     = @Unitary;                         %matrix type
p.modes      = 20;                               %matrix size 
p.part{1}    = 20;                               %bin size 
p.Max        = 40;                               %max counts
p.name       = sprintf('+P uniformsq, M=%d, Pa=%d',p.modes,p.part{1}); 
p.tr         = 0.5;                              %transmission
p.sqz        = ones(1,20);                       %uniform squeezing                               %decoherence factor
p.cutoff     = 1.e-6;                            %chi-2 cutoff
p.ensembles  = [10000,100,12];                   %repeats for errors
p.observe    = {@N1};                            %total number probability
p.compare    = {@N1lsc};                         %exact comparison
p.glabels    = {{'Counts m'}};                   %graphics x-label
p.olabels    = {'P(m)'};                         %graphics y-label
p.xk{1}      = {0:p.Max};                        %computed counts
p.diffplot   = {1};                              %plot the differences
[e1,d,p]    = xqsim(p);                          %run the simulation
xgraph(d,p);                                     %graph the data