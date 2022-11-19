function C = nm(a,p)
% C = NM(a,p) generates moments of photon number n. 
%  The moment list is defined by the vector O

O = p.O{p.k};                                    %list of moments
Ol = length(O);                                  %length of graph axis
sm  = (p.method-1.0)/2;                          %ordering correction
np2 = a(1:p.M,:).*a(p.M+1:2*p.M,:)-sm;           %normal ordered numbers
C = zeros(Ol,p.ensembles(1));                    %initialize correlation
for i = 1:Ol                                     %loop over length of graph
    C(i,:) = real(prod(np2(1:O(i),:),1));        %get O-th order moment
end                                              %end loop over length
C = mean(C,2);                                   %Average over ensemble
end                                              %end om function