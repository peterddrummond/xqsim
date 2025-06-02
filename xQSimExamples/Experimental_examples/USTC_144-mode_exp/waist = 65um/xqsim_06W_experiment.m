function e = xqsim_06W_experiment( )

p.name = '144-mode: Waist = 65um, power = 0.6W';  %Choose run name
p.matrix    = @expmatrix;                         %Choose matrix type
p.modes         = 144;                                %matrix size
p.N         = 50;                                 %squeezed input size
p.bits      = 32;                                 %bits for data extraction
p.sqz         = @expsqueeze;                        %squeezing parameter(s)                      
p.part{3} = [72, 72];                             %two-fold partition
p.ensembles = [5000,20,12];                      %long repeats for errors
% p.tr         = 1.0015;                             %transmission correction
% p.thermal       = 0.0260;                            %decoherence factor dc
p.t         = 1.00;                             %transmission correction
p.eps       = 0;                            %decoherence factor dc
p.tmss = 1;
% [Count]  = qcountsle(p);                        %data extraction
% p.counts = Count;                               %experimental coutns
p.counts    = 45782667;                           %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.observe   = {@K,@K1,@Kn};
p.compare   = {@expk,@expk1,@expk2};
p.logs{2}   = {0,1};
p.logs{3}   = {0,0,1};
p.diffplot  = {2,2,2};
p.glabels   = {{'Mode m'},{'Clicks m'},{'Clicks m_1','Clicks m_2'}};
p.olabels   = {'<G^{(1)}(m)>','G^{(M)}(m)','G^{(M)}(m)'};
p.xk{2}     = {0:p.modes};
p.xk{3}     = {0:72,0:72};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);
