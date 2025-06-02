function C = Y(a,p)
% C = Y(a,p) generates the y quadrature 

if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  C = -1i*(a(1:p.modes,:)- a(p.modes+1:2*p.modes,:));%y amplitude
else                                             % else not +P                          
  C = -1i*(a(1:p.modes,:)- a(1:p.modes,:));       %y amplitude
end                                              %end function