function C = km3(a,p)
% C = km3(a,p); generates ALL third-order click correlation moments for
% sequential modes where j<k<h. 


% INITIALIZE PROJECTORS

if p.method ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end

n = a(1:p.M,:).*a(p.M+1:2*p.M,:);                %normal ordered numbers
k  = 1-exp(-n);                                  %click projector

% COMPUTE CORRELATIONS

O = nchoosek(p.M,p.CO{p.k});
C = zeros(O,p.ensembles(1));
E = zeros(O,p.ensembles(1));
s=1;

for j = 1:p.M-2
    for h = j:p.M-2
        for i = h:p.M-2
            E(s,:) = real(k(j,:).*k(h+1,:).*k(i+2,:));
            s = s+1;
        end 
    end
end 

C = C + E;
C = mean(C,2);
end 