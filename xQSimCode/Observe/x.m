function C = x(a,p)
% C = X(a,p) generates the x quadrature for qsim

C = (a(1:p.M,:)+ a(p.M+1:2*p.M,:));           %y amplitude
C = mean(C,2);                                   %Average over ensemble
end                                              %end x function