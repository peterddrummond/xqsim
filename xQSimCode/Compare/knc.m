function C = knc(p)
% C = KNC(p); generates comparison click probabilities 
% for an n-fold partition. 
n = (sinh(p.r').*p.t').^2;
m = ((1-p.eps').*cosh(p.r').*sinh(p.r')).*(p.t'.^2);
k = 1-sqrt(1./((n+1).^2-m.^2));
if iscell(p.Part{p.k})                           %check if cell partition 
    C = xqchooseftnc(p.Part{p.k},k);              %get click probabilities
else                                             %not cell partition 
    C = xqchooseftn(p.Part{p.k},k);               %get click probabilities
end                                              %check if cell partition 
end 
