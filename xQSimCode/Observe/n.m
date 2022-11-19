function C = n(a,p)
% C = N(a,p); generates mean numbers of photons per channel for qsim

sm  = (p.method-1.0)/2;                     %ordering correction factor
C = a(1:p.M,:).*a(p.M+1:2*p.M,:)-sm;        %s-ordered numbers
C = mean(real(C),2);                        %take vector parallel mean
end                                         %End np function