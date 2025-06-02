function Kp = K1(a,varargin)
% Kp = K1test(a,p); generates click probabilities in a single partition

switch nargin
  case 4
    modes = varargin{1};
    xk   = varargin{2};
    p    = varargin{3};
  case 2
    p    = varargin{1};
    modes = p.part{p.nobserve}{1};
    xk   = p.xk{p.nobserve}{1};
  otherwise
    error('K1 requires 2 or 4 inputs, not %d', nargin);
end
if p.phase ~= 1
    error('xqsim supports click probabilities only for +P method = 1');
end
nmodes = modes(end)-modes(1)+1;
if isempty(xk)
    xk = 0:length(modes);
end
N0  = xk(1);
Nx  = xk(end);
N   = Nx-N0; ens = size(a,2);
iktheta = 2*pi*1i*(0:N)'/(N+1);
np2 = a(modes,:).*a(p.modes+modes,:);  %normal ordered numbers
cp  = exp(-np2);
exk = exp(-iktheta);
Kp  = ones(N+1,ens);
cp1 = 1-cp;
for i = 1:nmodes
  Kp = Kp.*(cp(i,:)+cp1(i,:).*exk);
end
Kp = Kp.*exp(iktheta*N0);                       %shift Fourier 
Kp = ifft(Kp);                                   %inverse Fourier transform
end                                              %End K1 function