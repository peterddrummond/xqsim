function C = Delv(a,p)
% C = DELV(a,p) generates the v quadrature standard deviation  for qsim
% follows conventions from the entanglement sction

fact =1/sqrt(p.m-1);
if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  y = (a(1:p.modes,:) - a(p.modes+1:2*p.modes,:))/1i;
else 
  y = (a(1:p.modes,:) - conj(a(1:p.modes,:)))/1i;% else not +P 
end                                              % end if positive P
v = y(1,:) + sum(y(2:p.m,:),1)*fact;
C = sqrt(real(mean(v.^2,2)) -2*(p.phase-2));
end                                              %end delv function