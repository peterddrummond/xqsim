function [e] = xqsim_entangle_short( )

%SET INPUT DATA PARAMETERS
p.name       = 'Multipartite-entanglement P200 ';     %Choose the run name
p.matrix     = @Beamsplitter;                    %Choose the matrix type
p.N          = 2;                                %squeezed input size N
p.Mvect          = 2:20:102;                              %matrix size M
p.ensembles  = [100,10,12];                     %repeats for errors
%p.ensembles  = [10,2,1];                     %repeats for errors
p.r          = [3,-3,0,0,0,0,0,0,0,0,0,0,0];     %squeezing parameter(s)
p.method     = 1;                                %method - Wigner
p.glabels    = {{'Modes m'},{'Modes m'},{'Modes m'},{'Modes m'}};
p.olabels    = {'\Delta u','\Delta v','(\Delta u)^2+(\Delta v)^2',...
               '\Delta u\Delta v '};
p.observe    = {@delu2,@delv2,@delsum,@delprod2};
p.compare{1} = @(p) (2*exp(-2*p.r(1)));
p.compare{2} = @(p) (2*exp(-2*p.r(1)));
p.compare{3} = @(p) 4*exp(-2*p.r(1));
p.compare{4} = @(p) 4*exp(-4*p.r(1));
p.xk{1}      = {p.Mvect};
p.xk{2}      = {p.Mvect};
p.xk{3}      = {p.Mvect};
p.xk{4}      = {p.Mvect};
j=1;
for m1 = p.Mvect                             %loop over matrix size
    p.M = m1;                                %set the matrix size
    [e,d,cp]  = xqsim(p);                       %repeat the simulation 
      for k = 1:4
        dat{1}{k}(1,j,:) = d{1}{k}(1,1,:);    %Store the output for method 
      end
      j=j+1;
end 
xgraph(dat,cp);