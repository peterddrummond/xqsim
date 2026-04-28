function e = xqsim_experiment_permute( )

%100-mode experiment with random permutation of the experimental data. 

%SET INPUT DATA PARAMETERS

p.name      = 'USTC 100-mode experiment - Random permutations';%Choose the run name
p.matrix    = @expmatrix;                        %Choose the matrix type
p.modes     = 100;                               %matrix size
p.sqz       = @expsqueeze;                       %squeezing parameter(s)
p.part{3}   = [50, 50];                          %two-fold partition
p.ensembles = [5000,20,12];                     %long repeats for errors
%p.ensembles = [1000,10,12];                       %repeats for errors
p.tr        = 1.0305;                            %transmission correction
p.thermal   = 0.0811;                             %decoherence factor dc
p.tmss = 1;                                     %Two-mode squeezing switch
p.counts    = 51392341;                          %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.permute   = xpermutation;
p.observe   = {@K,@K1,@Kn};
p.compare   = {@expk,@expk1,@expk2};
p.logs{2}   = {0,1};
p.logs{3}   = {0,0,1};
p.diffplot  = {2,2,2};
p.glabels   = {{'j'},{'m'},{'Clicks m_1','Clicks m_2'}};
p.olabels   = {'<\pi_j(1)>','G(m)','G_2(m)'};
p.xk{2}     = {0:p.modes};
p.xk{3}     = {0:p.modes/2,0:p.modes/2};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);