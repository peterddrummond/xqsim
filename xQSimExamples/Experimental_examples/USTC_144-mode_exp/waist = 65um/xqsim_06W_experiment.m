function e = xqsim_06W_experiment( )

p.name = '144-mode: Waist = 65um, power = 0.6W';  %Choose run name
p.matrix    = @expmatrix;                         %Choose matrix type
p.M         = 144;                                %matrix size
p.N         = 50;                                 %squeezed input size
p.bits      = 32;                                 %bits for data extraction
p.r         = @expsqueeze;                        %squeezing parameter(s)                      
p.Part{3} = [72, 72];                             %two-fold partition
%p.ensembles = [5000,20,12];                      %long repeats for errors
p.ensembles = [1000,10,12];                       %repeats for errors
p.t         = 0.9966;                             %transmission correction
p.eps       = 0.02889;                            %decoherence factor dc
% [Count]  = qcountsle(p);                        %data extraction
% p.counts = Count;                               %experimental coutns
p.counts    = 18523709;                           %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@k,@k1,@kn};
p.compare   = {@expk,@expk1,@expk2};
p.logs{2}   = [0,1];
p.logs{3}   = [0,0,1];
p.diffplot  = {2,2,2};
p.glabels   = {{'Mode m'},{'Clicks m'},{'Clicks m_1','Clicks m_2'}};
p.olabels   = {'<G^{(1)}(m)>','G^{(M)}(m)','G^{(M)}(m)'};
p.xk{2}     = {0:p.M};
p.xk{3}     = {0:72,0:72};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);
