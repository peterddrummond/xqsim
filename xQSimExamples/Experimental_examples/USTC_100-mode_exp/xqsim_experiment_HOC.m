function e = xqsim_experiment_HOC( )

% HOC = High-order correlations
% This algorithm is a test for higher order correlations for 100-mode data. 

p.name      = 'USTC 100-mode experiment - HOC';  %Choose the run name
p.matrix    = @expmatrix;                        %Choose the matrix type
p.modes     = 100;                               %matrix size
p.sqz       = @expsqueeze;                       %squeezing parameter(s)
p.CO{2}     = 2;
p.CO{3}     = 2;
p.ensembles = [600,200,12];
%p.ensembles = [6000,200,12];
p.tr        = 1.0305;                            %transmission correction
p.thermal   = 0.0811;                             %decoherence factor dc
p.tmss = 1;                                     %Two-mode squeezing switch
p.counts    = 51392341;                          %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@K,@Kmsub,@Km2};
p.compare   = {@expk,@expk2ms,@expk2m};
p.diffplot  = {2,2,2};
p.glabels   = {{'j'},{'j,k'},{'j<k'}};
p.olabels   = {'<\pi_j(1)>','<\pi_j(1)\pi_k(1)>','<\pi_j(1)\pi_k(1)>'};
p.xk{2}     = {0:p.modes-1};
p.xk{3}     = {0:nchoosek(p.modes,p.CO{3})};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);