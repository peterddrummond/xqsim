function C = Delsum(a,p)
% C = Delv(a,p) generates the sum of u,v variances for qsim
% follows conventions from the entanglement sction

C = Delv2(a,p)+Delu2(a,p);
end                                              %end du function