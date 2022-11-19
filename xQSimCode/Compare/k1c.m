function C = k1c(p)
% C = K1C(p); generates comparison click probabilities in one partition,
% assumed to be of length p.N

n = (sinh(p.r').*p.t').^2;
m = ((1-p.eps').*cosh(p.r').*sinh(p.r')).*(p.t').^2;
k = 1-sqrt(1./((n+1).^2-m.^2));
C = xqchooseft1(p.N, k);
C = C';
C(p.N+1:p.M+1,1) = 0; 
end 
