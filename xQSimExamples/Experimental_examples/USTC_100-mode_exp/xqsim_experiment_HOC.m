function e = xqsim_experiment_HOC( )

% HOC = High-order correlations
% This algorithm is a test for higher order correlations for 100-mode data. 

p.name      = 'Exp: ens=1.2e6, t=1.0235, e=0.0932'%Choose the run name
p.matrix    = @expmatrix;                        %Choose the matrix type
p.M         = 100;                               %matrix size
p.N         = 50;                                %squeezed input size
p.r         = @expsqueeze;                       %squeezing parameter(s)
p.CO{2}     = 2;
p.CO{3}     = 3;
p.CO{4}     = 4;
p.CO{5}     = 2;
p.ensembles = [6000,200,12];
p.t         = 1.0235;                            %transmission correction
p.eps       = 0.0932;                             %decoherence factor dc
p.counts    = 51392341;                          %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@k,@kmsub,@kmsub,@kmsub,@km2};
p.compare   = {@expk,@expk2ms,@expk3ms,@expk4ms,@expk2m};
p.diffplot  = {2,2,2,2,2};
p.glabels   = {{'j'},{'j,k'},{'j,k,h'},{'j,k,h,i'},{'j<k'}};
p.olabels   = {'<\pi_j(1)>','G^{(2)}','G^{(3)}','G^{(4)}','<\pi_j(1)\pi_k(1)>'};
p.xk{2}     = {0:p.M-1};
p.xk{3}     = {0:p.M-2};
p.xk{4}     = {0:p.M-3};
p.xk{5}     = {0:4950};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);