function C = nmc(p)
% C = NMC(p); generates comparison O-th order moments of photon number

O = p.O{p.k};                                    %list of moments
Ol = length(O);
n = (sinh(p.r').*p.t').^2;
for i = 1:Ol
    C(i,1) = prod(n(1:O(i)));                    %Get O-th order moment
end
end 