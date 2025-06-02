function C = Delu(a,p)
% C = delu(a,p) generates the u quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.m-1);

if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  x = a(1:p.m,:) + a(p.m+1:2*p.m,:);
else 
  x = (a(1:p.modes,:) + conj(a(1:p.modes,:)));% else not +P 
end                                              % end if positive P
u = x(1,:) - sum(x(2:p.m,:),1)*fact;
C = sqrt(real(mean(u.^2,2)) -2*(p.phase-2));
end                                              %end du function