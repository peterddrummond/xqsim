function e = xqsim_141W_experiment_HOC( )

p.name = '144-mode: Waist = 125um, power = 1.41W'; %Choose run name
p.matrix    = @expmatrix;                          %Choose matrix type
p.M         = 144;                                 %matrix size
p.N         = 50;                                  %squeezed input size
p.r         = @expsqueeze;                         %squeezing parameter(s)
p.CO{2}     = 2;                                   %correlation order
p.CO{3}     = 3;                                   %correlation order
p.CO{4}     = 4;                                   %correlation order
p.ensembles = [6000,200,12];                        %long repeats for errors   
p.t         = 0.9972;                              %transmission correction
p.eps       = 0.0354;                              %decoherence factor dc
p.counts    = 46425243;                            %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.esample  = {-2,-2};
p.observe   = {@k,@kmsub,@kmsub,@kmsub};
p.compare   = {@expk,@expk2ms,@expk3ms,@expk4ms};
p.diffplot  = {2,2,2,2};
p.glabels   = {{'j'},{'j,k'},{'j,k,h'},{'j,k,h,i'}};
p.olabels   = {'<\pi_j(1)>','G^{(2)}','G^{(3)}','G^{(4)}'};
p.xk{2}     = {0:p.M-1};
p.xk{3}     = {0:p.M-2};
p.xk{4}     = {0:p.M-3};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);
