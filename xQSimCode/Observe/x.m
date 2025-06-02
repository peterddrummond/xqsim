function C = X(a,p)
% C = X(a,p) generates the x quadrature 

if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  C = (a(1:p.modes,:)+ a(p.modes+1:2*p.modes,:));%x amplitude
else                                             % else not +P 
  C = a(1:p.modes,:)+ conj(a(1:p.modes,:));      %x amplitude
end                                              %end x function