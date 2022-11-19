function [U] = Unitary(p)
%[U] = UNITARY(M)
%Calculates a random complex unitary matrix, size M*M

if isnumeric(p)  && p <= 0
    U =  'complex unitary';
    return;
end
U = (randn(p.M) + 1i*randn(p.M))/(sqrt(2));
[U,R] = qr(U);
D = diag(R);
ph = D./abs(D);
U = U*diag(ph);
end