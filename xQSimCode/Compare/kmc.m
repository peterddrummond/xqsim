function C = kmc(p)
% C = KMC(p); generates comparison O-th order click probability 

O = p.O{p.k};                                    %list of moments
Ol = length(O);                                  %length of graph axis 
n = (sinh(p.r').*p.t').^2;
m = ((1-p.eps').*cosh(p.r').*sinh(p.r')).*(p.t').^2;
cp = 1-sqrt(1./((n+1).^2-m.^2));
C = zeros(Ol,1);                                 %initialize correlation
for i = 1:Ol
    C(i,1) = prod(cp(1:O(i)));                   %Get O-th order click prob
end
end 
