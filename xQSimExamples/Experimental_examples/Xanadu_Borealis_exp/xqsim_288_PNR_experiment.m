function e = xqsim_288_PNR_experiment()

p.name = '288-mode GBS with PNR'; 
p.matrix    = @expmatrix;                        
p.modes         = 288;                               
p.N         = 288;
p.max       = 300;                    %Upper bound on number of recordered counts
p.sqz         = @expsqueeze;  
p.ensembles = [5000,20,12];
p.tr         = 1.00;                              %transmission factor
p.thermal       = 0;                              %decoherence factor dc
% p.tr         = 0.9848;                              %transmission factor
% p.thermal       = 0.0547;                              %decoherence factor dc
% [Count]  = xqcountspnr(p);
% p.counts = Count;
p.counts    = 1113000;                          %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@N1};
p.compare   = {@expn1};
p.logs{1}   = {0,1};
p.diffplot  = {2};
p.glabels   = {{'m'}};
p.olabels   = {'G^{(M)}(m)'};
p.xk{1}     = {0:p.max};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);
