function e = xqsim_03W_experiment( )

p.name = '144-mode: Waist = 65um, power = 0.3W';  %Choose run name
p.matrix    = @expmatrix;                         %Choose matrix type
p.M         = 144;                                %matrix size
p.N         = 50;                                 %squeezed input size
p.bits      = 32;                                 %bits for data extraction
p.r         = @expsqueeze;                        %squeezing parameter(s)                      
p.Part{3} = [72, 72];                             %two-fold partition
%p.ensembles = [5000,20,12];                      %long repeats for errors
p.ensembles = [1000,10,12];                       %repeats for errors
p.t         = 0.9972;                             %transmission correction
p.eps       = 0.0202;                             %decoherence factor dc
% [Count]  = qcountsle(p);                        %data extraction
% p.counts = Count;                               %experimental counts
p.counts    = 46535311;                           %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@k,@k1,@kn};
p.compare   = {@expk,@expk1,@expk2};
p.logs{2}   = [0,1];
p.logs{3}   = [0,0,1];
p.diffplot  = {2,2,2};
p.glabels   = {{'j'},{'m'},{'Clicks m_1','Clicks m_2'}};
p.olabels   = {'<\pi_j(1)>','G^{(M)}(m)','G^{(M)}(m)'};
p.xk{2}     = {0:p.M};
p.xk{3}     = {0:72,0:72};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);
