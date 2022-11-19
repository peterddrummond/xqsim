function C = delu2(a,p)
% C = delu(a,p) generates the u quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.M-1);
x = a(1:p.M,:) + a(p.M+1:2*p.M,:);
u = x(1,:) - sum(x(2:p.M,:),1)*fact;
C = real(mean(u.^2,2)) -2*(p.method-2);
end                                              %end du function