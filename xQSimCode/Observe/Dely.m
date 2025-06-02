function C = Dely(a,p)
% C = Dely(a,p) generates the y quadrature standard deviation

C = 2 - p.phase;                                 %ordering correction
if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  C = C - (a(1:p.m,:) - a(p.m+1:2*p.m,:)).^2;    %Y amplitude squared
else                                             % else not +P 
  C = C - (a(1:p.m,:) - conj(a(1:p.m,:))).^2;    %Y amplitude squared
end                                              % end if positive P
C = sqrt(real(mean(C,2)));                       %Average over ensemble
end                                              %end x function