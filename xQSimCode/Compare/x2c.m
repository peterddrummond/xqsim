function C = X2c(p)
% C = X2C(a,p); generates comparison quadrature moments per channel
% assumes output unchanged from input. 
% eg, identity matrix OR uniform thermal+unitary

n = (sinh(p.sqz').*p.tr').^2;
m = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz')).*(p.tr').^2;
C = 2*(n+m)+1;
C = reshape(C,[1,1,numel(C)]);
end                                              %end x2c function