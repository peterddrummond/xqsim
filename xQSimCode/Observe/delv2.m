function C = Delv2(a,p)
% C = DELV2(a,p) generates the v quadrature variance  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.modes-1);
if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  y = (a(1:p.modes,:) - a(p.modes+1:2*p.modes,:))/1i;
else 
  y = (a(1:p.modes,:) - conj(a(1:p.modes,:)))/1i;% else not +P 
end                                              % end if positive P
v = y(1,:) + sum(y(2:p.modes,:),1)*fact;
C = real(mean(v.^2,2)) -2*(p.phase-2);
end                                              %end du function