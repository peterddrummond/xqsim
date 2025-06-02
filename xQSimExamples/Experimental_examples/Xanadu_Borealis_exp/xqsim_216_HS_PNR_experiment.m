function e = xqsim_216_HS_PNR_experiment()

p.name = '216-mode GBS with PNR - High squeeze'; 
p.matrix    = @expmatrix;                        
p.modes         = 216;                               
p.N         = 216;
p.max       = 250;                    %Upper bound on number of recordered counts
p.sqz         = @expsqueeze;   
p.ensembles = [5000, 20, 12];
p.tr         = 1.00;                              %transmission factor
p.thermal       = 0;                              %decoherence factor dc
% p.tr         = 0.9941;                              %transmission factor
% p.thermal       = 0.0510;                              %decoherence factor dc
% [Count]  = xqcountspnr(p);
% p.counts = Count;
p.counts    = 1026000;                          %experimental counts
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
