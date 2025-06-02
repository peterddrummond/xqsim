function C = Knc(p)
% C = KNC(p); generates comparison click probabilities 
% for an n-fold partition. Uses limits defined by p.xk
n = (sinh(p.sqz').*p.tr').^2;
m = ((1-p.thermal').*cosh(p.sqz').*sinh(p.sqz')).*(p.tr').^2;
k = 1-sqrt(1./((n+1).^2-m.^2));
if iscell(p.part{p.noutput})                     %check if cell partition 
    C = xchooseftnc(p.part{p.noutput},k);        %get click probabilities
else                                             %not cell partition 
    C = xchooseftn(p.part{p.noutput},k);         %get click probabilities
end                                              %check if cell partition
in = [p.xk{p.noutput}{:},{0},{0},{0},{0},{0}];
C = C(1+in{1}(1):1+in{1}(end),1+in{2}(1):1+in{2}(end),...
    1+in{3}(1):1+in{3}(end),1+in{4}(1):1+in{4}(end),...
    1+in{5}(1):1+in{5}(end),1+in{6}(1):1+in{6}(end)); 
C = reshape(C,[1,1,size(C)]);
end 