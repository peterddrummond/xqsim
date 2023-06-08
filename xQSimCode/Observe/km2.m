function C = km2(a,p)
% C = km2(a,p); generates ALL second-order click correlation moments for
% sequential modes where j<k. 


% INITIALIZE PROJECTORS

if p.method ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end

n = a(1:p.M,:).*a(p.M+1:2*p.M,:);                %normal ordered numbers
k  = 1-exp(-n);                                  %click projector

% COMPUTE CORRELATIONS

O = nchoosek(p.M,p.CO{p.k});
C = zeros(O,p.ensembles(1));
A = cell(1,p.M-1);


for i=1:p.M-1
    for j=i:p.M-1
        A{i}(j-i+1,:) = real(k(i,:).*k(j+1,:));
    end 
end 
C = C + cat(1,A{:});
C = mean(C,2);
end 