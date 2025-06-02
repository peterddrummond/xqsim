function C = Nm(a,p)
% C = NM(a,p) generates correlations of photon number n.
% The moment index list is defined by the vector correl.
% Duplicate mode indices are not allowed, unless normal-ordered

ens = size(a,2);
O = p.correl{p.nobserve};                        %define correlations
Ol = length(O);                                  %length of graph axis
sm  = (p.phase-1.0)/2;                           %ordering correction
if p.phase == 1 && size(a,1) >= 2*p.modes        % if positive P
  N = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);  % normal-ordered numbers
else                                             % else not +P 
  N = a(1:p.modes,:).*conj(a(1:p.modes,:))-sm;   % add ordering correction
end                                              % end if positive P
C = zeros(Ol,ens);                               %initialize correlation
for i = 1:Ol                                     %loop over length of graph
    C(i,:) = real(prod(N(1:O(i),:),1));          %get O-th order moment
end                                              %end loop over length
end                                              %end om function