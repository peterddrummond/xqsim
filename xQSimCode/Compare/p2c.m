function C = p2c(p)
% C = P2C(a,p); generates comparison P quadrature moments per channel
% assumes output unchanged from input. 
% eg, identity matrix or uniform thermal+unitary

n = (sinh(p.r').*p.t').^2;
m = ((1-p.eps').*cosh(p.r').*sinh(p.r')).*(p.t').^2;
C = 2*(n-m)+1;
end                                              %end y2c function