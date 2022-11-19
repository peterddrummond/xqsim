function C = k(a,p)
% C = k(a,p); generates mean numbers of clicks per channel for qsim

if p.method ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end
n = a(1:p.M,:).*a(p.M+1:2*p.M,:);                %normal ordered numbers
C  = 1-exp(-n);                                  %click projector
C = mean(C,2);                                   %Average over ensemble
end                                              %End nc function