function e1 = xy_quad_psqueeze( )
%Script demonstrates binranges by computing the bivariate x,y-quadrature
%probability for positive-P amplitude inputs. These probabilities are
%compared to exact densisites pure squeezed state inputs. Valid
%for xqsim3.0.

p.matrix         = @Identity;                      %matrix type
p.modes          = 1;                              %matrix size m
p.name           = sprintf('+P pure squeezed input state, M=%d',p.modes); 
p.tr             = 1;                              %transmission
p.sqz            = ones(1,p.modes);                %uniform squeezing
p.thermal        = zeros(1,p.modes);               %decoherence factor 
p.binranges{1}   =  {-5:0.25:5,-5:0.25:5};         %binning ranges
p.ensembles      = [10000,100,12];                 %repeats for errors
p.cutoff         = 1.e-5;                          %output cutoff
p.observe        = {@Gauss2sqz};                   %observe code
p.compare        = {@Gauss2sqzc};                  %compare code
p.glabels        = {{'x','y'}};                    %horizontal axis labels
p.olabels        = {'P(x,y)'};                     %vertical axis labels
p.diffplot       = {2};                            %plot difference
p.xk{1}          = {-5:0.25:5,-5:0.25:5};          %probability axis range
[e1,d,cp]        = xqsim(p);
xgraph(d,cp);

% --------------- joint x,y-quadrature observe code ---------------

function C = Gauss2sqz(a,p)
    if ~isequal(p.thermal, zeros(1,p.modes))
        error('pure squeezed inputs required');
    end
    x = a(1:p.modes,:) + a(p.modes+1:2*p.modes,:);
    y = (a(1:p.modes,:) - a(p.modes+1:2*p.modes,:));
    C = [x;y];
end 

% ----- Exact pure squeezed state bivariate binned prob density ------

function C = Gauss2sqzc(p)
    if ~isequal(p.thermal, zeros(1,p.modes))
        error('pure squeezed inputs required');
    end
    np = (sinh(p.sqz')).^2;
    mp = cosh(p.sqz').*sinh(p.sqz');
    sigx = 2*(np + mp);
    sigy = 2*(mp - np);
    [x,y] = ndgrid(p.oc{1}{1},p.oc{1}{2});
    nm = 1./sqrt(4*pi^2*sigx*sigy);
    C = nm.*exp(-x.^2/(2*sigx));
    C = C.*exp(-y.^2./(2*sigy));
end 

end 
    