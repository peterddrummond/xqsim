function e1 = x_quad_1D_marginal( )
%Script demonstrates binranges by computing the 1D marginal x-quadrature
%probability for positive-P amplitude inputs. These probabilities are
%compared to exact marginals for thermal and pure squeezed states. Valid
%for xqsim3.0.

p.matrix       = @Identity;                        %matrix type
p.modes        = 1;                                %matrix size m
p.name         = sprintf('+P thermal input state, M=%d',p.modes); 
p.tr           = 1;                                %transmission
p.sqz          = ones(1,p.modes);                  %uniform squeezing
p.thermal      = ones(1,p.modes);                  %decoherence factor 
p.binranges{1} =  {-5:0.25:5};                     %binning ranges
p.ensembles    = [1000,100,12];                    %repeats for errors
p.cutoff       = 1.e-7;                            %output cutoff
p.observe      = {@Gauss1};                        %observe code
p.compare      = {@Gauss1c};                       %compare code
p.glabels      = {{'x'}};                          %Horizontal axis label
p.olabels      = {'P(x)'};                         %vertical axis label
p.diffplot     = {2};                              %plot difference
p.xk{1}        = {-5:0.25:5};                      %probability axis range

p1 = p;                                            %secondary data
p1.thermal = zeros(1,p.modes);                     %decoherence factor
p1.name= sprintf('+P pure squeezed input state, M=%d',p.modes); 

[e1,d,cp]    = xqsim(p,p1);
xgraph(d,cp);


% --------------- x-quadrature observe code ---------------

function C = Gauss1(a,p)
    %Test is for thermal or squeezed states
    x = a(1:p.modes,:) + a(p.modes+1:2*p.modes,:);%x-quadrature
    C = x;
end 

% ----------- Exact x-quad marginal probability density --------------

function C = Gauss1c(p)
    %Thermal/Squeezed state distribution test 
    np = (sinh(p.sqz')).^2;
    mp = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz'));
    sigx = 2*(np + mp);
    x = p.oc{1}{1};
    nm = 1./sqrt(4*pi*(np+mp));
    nm = nm(1);
    C = nm.*exp(-x.^2./(2*sigx(1))); 
end

end 
    