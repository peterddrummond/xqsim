function C = delv2(a,p)
% C = delv2(a,p) generates the v quadrature variance  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.M-1);
y = (a(1:p.M,:) - a(p.M+1:2*p.M,:))/1i;
v = y(1,:) + sum(y(2:p.M,:),1)*fact;
C = real(mean(v.^2,2)) -2*(p.method-2);
end                                              %end du function