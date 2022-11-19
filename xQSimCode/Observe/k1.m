function C = k1(a,p)
% C = k1(a,p); generates click probabilities in a single partition

if p.method ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end
np2 = a(1:p.M,:).*a(p.M+1:2*p.M,:);         %normal ordered numbers
cp  = 1-exp(-np2);                          %click projector
C = xqchooseft1(p.M,cp);
sz = size(C);
C = reshape(C, sz(2), sz(1));
end                                         %End cprob1 function