function C = Delu2(a,p)
% C = DELU2(a,p) generates the u quadrature variance  for qsim
% follows conventions from the entanglement section
fact =1/sqrt(p.modes-1);                     
if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  x = a(1:p.modes,:) + a(p.modes+1:2*p.modes,:);
else                                             % else not +P 
  x = a(1:p.modes,:) + conj(a(1:p.modes,:));
end                                              % end if positive P
u = x(1,:) - sum(x(2:p.modes,:),1)*fact;
C = real(mean(u.^2,2)) -2*(p.phase-2);
end                                              %end du function