function C = p(a,p)
% C = P(a,p) generates the p quadrature for qsim

C = -i*(a(1:p.M,:)+ a(p.M+1:2*p.M,:));           %y amplitude
C = mean(C,2);                                   %Average over ensemble
end                                              %end x function