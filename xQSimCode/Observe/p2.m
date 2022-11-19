function C = p2(a,p)
% C = P2(a,p) generates the x quadrature second moment for qsim

C = 2 - p.method;                                %ordering correction
C = C - (a(1:p.M,:)- a(p.M+1:2*p.M,:)).^2;       %x amplitude squared
C = mean(C,2);                                   %Average over ensemble
end                                              %end x function