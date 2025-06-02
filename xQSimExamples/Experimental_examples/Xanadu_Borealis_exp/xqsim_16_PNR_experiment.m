function e = xqsim_16_PNR_experiment()

p.name = '16-mode GBS with PNR'; 
p.matrix    = @expmatrix;                        
p.modes         = 16;                               
p.N         = 16;
p.max       = 50;                    %Upper bound on number of recordered counts
p.sqz         = @expsqueeze; 
p.ensembles = [5000, 20, 12];
p.tr         = 1.00;                              %transmission factor
p.thermal       = 0;                              %decoherence factor dc
% p.tr         = 0.9841;                           %transmission factor
% p.thermal       = 0.0648;                           %decoherence factor dc
% [Count]  = xqcountspnr(p);
% p.counts = Count;
p.counts    = 84096000;                          %experimental counts
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