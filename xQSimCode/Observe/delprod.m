function C = Delprod(a,p)
% C = Delprod(a,p) generates the product of u,v standard deviations
% follows conventions from the entanglement sction

C = Delv(a,p).*Delu(a,p);
end                                              %end du function