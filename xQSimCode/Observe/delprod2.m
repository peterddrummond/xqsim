function C = Delprod2(a,p)
% C = Delprod(a,p) generates the product of u,v variances
% follows conventions from the entanglement sction

C = Delv2(a,p).*Delu2(a,p);
end                                              %end du function