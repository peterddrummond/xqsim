function C = kn(a,p)
% C = KN(a,p); generates click probabilities for an n-fold partition
% The partition is defined by the quantity p.Part{p.k}
% This is EITHER  a cell array of index vectors
% OR, a vector of partition lengths

if p.method ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end
n = a(1:p.M,:).*a(p.M+1:2*p.M,:);                %normal ordered numbers
k  = 1-exp(-n);                                  %click projector
if iscell(p.Part{p.k})                           %check if cell partition 
    C = xqchooseftnc(p.Part{p.k},k);              %get click probabilities
else                                             %not cell partition 
    C = xqchooseftn(p.Part{p.k},k);               %get click probabilities
end                                              %check if cell partition 
end                                              %end kn