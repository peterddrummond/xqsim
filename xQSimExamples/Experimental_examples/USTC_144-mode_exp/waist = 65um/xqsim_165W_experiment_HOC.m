function e = xqsim_165W_experiment_HOC( )

p.name = '144-mode: Waist = 65um, power = 1.65W'; %Choose run name
p.matrix    = @expmatrix;                         %Choose matrix type
p.M         = 144;                                %matrix size
p.N         = 50;                                 %squeezed input size
p.r         = @expsqueeze;                        %squeezing parameter(s) 
p.CO{2}     = 2;                                  %correlation order
p.CO{3}     = 3;                                  %correlation order
p.CO{4}     = 2;                                  %correlation order
p.CO{5}     = 3;                                  %correlation order
p.CO{6}     = 4;                                  %correlation order
p.ensembles = [6000,200,12];                      %long repeats for errors
p.t         = 1.01099;                            %transmission correction
p.eps       = 0.04282;                            %decoherence factor dc
p.counts    = 40183178;                           %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@k,@km2,@km3,@kmsub,@kmsub,@kmsub};
p.compare   = {@expk,@expk2m,@expk3m,@expk2ms,@expk3ms,@expk4ms};
p.diffplot  = {2,2,2,2,2,2};
p.glabels   = {{'j'},{'j<k'},{'j<k<h'},{'j,k'},{'j,k,h'},{'j,k,h,i'}};
p.olabels   = {'<\pi_j(1)>','<\pi_j(1)\pi_k(1)>','<\pi_j(1)\pi_k(1)\pi_h(1)>>','G^{(2)}','G^{(3)}','G^{(4)}'};
p.xk{2}     = {0:10296};
p.xk{3}     = {0:487344};
p.xk{4}     = {0:p.M-1};
p.xk{5}     = {0:p.M-2};
p.xk{6}     = {0:p.M-3};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);
