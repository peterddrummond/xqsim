function C = delprod2(a,p)
% C = delprod(a,p) generates the product of u,v variances
% follows conventions from the entanglement sction

C = delv2(a,p).*delu2(a,p);
end                                              %end du function