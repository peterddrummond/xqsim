function e = xqsim_experiment_permute( )

%100-mode experiment with random permutation of the experimental data. 

%SET INPUT DATA PARAMETERS

p.name      = 'Exp: ens=1.2e6, t=1.0235, e=0.0932'%Choose the run name
p.matrix    = @expmatrix;                        %Choose the matrix type
p.M         = 100;                               %matrix size
p.N         = 50;                                %squeezed input size
p.r         = @expsqueeze;                       %squeezing parameter(s)
p.Part{3}   = [50, 50];                          %two-fold partition
p.Part{4}   = [25,25,25,25];
p.ensembles = [250,400,12];
%p.ensembles = [5000,20,12];                     %long repeats for errors
%p.ensembles = [1000,10,12];                       %repeats for errors
p.t         = 1.0235;                            %transmission correction
p.eps       = 0.0932;                             %decoherence factor dc
p.counts    = 51392341;                          %experimental counts
p.mincount  = 10;
p.cutoff    = 1e-7;
p.permute   = xqpermutation;
p.observe   = {@k,@k1,@kn,@kn};
p.compare   = {@expk,@expk1,@expk2,@expk4};
p.logs{2}   = [0,1];
p.logs{3}   = [0,0,1];
p.logs{4}   = [0,0,0,1];
p.diffplot  = {2,2,2,2};
p.glabels   = {{'j'},{'m'},{'Clicks m_1','Clicks m_2'},...
              {'Clicks m_1','Clicks m_2','Clicks m_3','Clicks m_4'}};
p.olabels   = {'<\pi_j(1)>','G(m)','G(m)','G(m)'};
p.xk{2}     = {0:p.M};
p.xk{3}     = {0:50,0:50};
p.xk{4}     = {0:25,0:25,0:25,0:25};
[e,d,cp]    = xqsim(p);
xgraph(d,cp);