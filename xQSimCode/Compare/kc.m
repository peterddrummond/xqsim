function C = kc(s)
% C = KC(a,p); generates comparison numbers of clicks per channel
% assumes output unchanged from input. 
% eg, identity matrix or uniform thermal+unitary

n = (sinh(s.r').*s.t').^2;
m = ((1-s.eps').*cosh(s.r').*sinh(s.r')).*(s.t').^2;
C=1-sqrt(1./((n+1).^2-m.^2));
end                                         %End ncc function