function [U] = Beamsplitter(p)
%[U] = BEAMSPLITTER(p)
%Calculates a generalized beamsplitter matrix, size p.m x p.m

if isnumeric(p)  && p <= 0
    U =  'beamsplitter';
    return;
end
M=p.M;
r = zeros(1,M);
t = zeros(1,M);
r1 = zeros(1,M+1);
r1(1) = -1;
r(1)=1/sqrt(2);
U = zeros(M,M);
for k =2:M
	r(k)=1/sqrt(M-k+1);
end
for k =1:M
	t(k)=sqrt(1-r(k)^2);
	r1(k+1)=r(k);
	if k<M
		U(k,k+1)=t(k);
	end
    for j =1:k
    		U(k,j) = -r(k)*prod(t(j:k-1))*r1(j);
    end
end
%test_for_unitarity = sum(sum((U*U'-eye(M))^2))
end