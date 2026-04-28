function e = xqsim_experiment_bestfit( )

%SET INPUT DATA PARAMETERS

p.name      = 'USTC 100-mode experiment';       %Choose the run name
p.matrix    = @expmatrix;                       %Choose the matrix type
p.modes     = 100;                              %matrix size
p.sqz       = @expsqueeze;                      %squeezing parameter(s)
p.part{3}   = [50, 50];                         %two-fold partition
p.ensembles = [500,20,12];                       %repeats for errors
p.ensembles = [5000,20,12];
p.tr        = 1.0305;                           %transmission correction
p.thermal   = 0.0811;                           %decoherence factor dc
p.tmss = 1;                                     %Two-mode squeezing switch
p.counts    = 51392341;                         %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@K,@K1,@Kn};
p.compare   = {@expk,@expk1,@expk2};
p.logs{2}   = {0,1};
p.logs{3}   = {0,0,1};
p.diffplot  = {2,2,2};
p.glabels   = {{'Mode m'},{'Clicks m'},{'Clicks m_1','Clicks m_2'}};
p.olabels   = {'G^{(1)}(m)>','G^{(M)}(m)','G^{(M)}(m)'};
p.xk{2}     = {0:p.modes};
p.xk{3}     = {0:p.modes/2,0:p.modes/2};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);