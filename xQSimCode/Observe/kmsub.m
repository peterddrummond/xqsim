function C = kmsub(a,p)

% C = KMSUB(a,p); Generates all a subset of all possible combinations of
% output mode click correlations for moments of sequential modes. This is
% done for graphical simplificity. 
% Example, for a second order correlation, the first data point is a 
% product of clicks from modes 1, 2, the second data point is a product 
% of clicks from modes 2,3, third data point is a product of clicks from 
% modes 3,4 etc.

% INITIALIZE PROJECTORS

if p.method ~= 1                                 %check method is OK     
    error('qSIM supports click probabilities only for +P method = 1');
end                                              %end check method is OK
               
n = a(1:p.M,:).*a(p.M+1:2*p.M,:);                %get normal ordered number
k  = 1-exp(-n);

% Compute correlations

Ord = p.CO{p.k} - 1;
NM = p.M - Ord;

C = zeros(NM, p.ensembles(1));

for i = 1:NM
    C(i,:) = real(prod(k(i:i+Ord,:),1));
end 
C = mean(C,2);                                   %Average over ensemble
end