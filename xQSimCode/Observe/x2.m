function C = X2(a,p)
% C = X2(a,p) generates the x quadrature second moment for xqsim

C = 2 - p.phase;                                 %ordering correction
if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  C = C + (a(1:p.modes,:)+ a(p.modes +1:2*p.modes,:)).^2; % +P case
else                                             % else not +P
  C = C + (a(1:p.modes,:)+ conj(a(1:p.modes,:))).^2; 
end                                              % end if normal-ordered 
end                                              %end x function