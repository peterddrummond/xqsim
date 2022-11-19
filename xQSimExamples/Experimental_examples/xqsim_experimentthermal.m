function e = xqsim_experimentthermal( )

%SET INPUT DATA PARAMETERS

p.name      = 'Exp: 1.2e6, t=1.0061, e=0.138';   %Choose the run name
p.matrix    = @expmatrix;                        %Choose the matrix type
p.M         = 100;                               %matrix size
p.N         = 50;                                %squeezed input size
p.r         = @expsqueeze;                       %squeezing parameter(s)
p.Part{3}   = [50, 50];                          % Two-fold partition
p.ensembles = [1000,100,12];                      %repeats for errors
%p.ensembles = [100,20,1];                      %repeats for errors
p.t         = 1.0061;                           %transmission factor
p.eps       = 0.138;                              %decoherence factor dc
p.counts    = 51392341;                          %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.esample  = {-2,-2};
p.observe   = {@k,@k1,@kn};
p.compare   = {@expk,@expk1,@expk2};
p.logs{2}   = [0,1];
p.logs{3}   = [0,0,1];
p.diffplot  = {2,2,2};
p.glabels   = {{'Mode m'},{'Clicks m'},{'Clicks m_1','Clicks m_2'}};
p.olabels   = {'G^{(1)}(m)>','G^{(M)}(m)','G^{(M)}(m)'};
p.xk{2}     = {0:p.M};
p.xk{3}     = {0:50,0:50};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);