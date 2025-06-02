function a = xqgen(a,p)
% a  = xqgen(a,p) generates phase-space data for boson sampling
% Input is a complex matrix, a = alpi;beti and a parameter struct p
% First the input amplitude is scaled by the recycling factor p.re
% This is combined with a local squeezed state or thermal excitation
% Squeezed state parameter is p.sqz, decoherence parameter is p.thermal 
% p.matrix is the transmission, p.phase is the phase-space phase = 1,2,3 
% Transmission is modified by transmission correction factor p.tr 
% Output is a = alp;bet;ngen, including decoherence and transmission factors
% ngen is the generated photon number 
% First a matrix index gives the channel count, repeated three times
% Second a matrix index gives the index of the parallel ensemble
% *Licensed by Peter D. Drummond, (2025) - see License.txt, xQSim manual  

%GENERATE PHASE-SPACE SAMPLES FOR ANY ORDERING

ens = p.ensembles(1);                            %ensemble number
sm  = (p.phase-1)/2;                             %ordering parameter
r   = p.sqz';                                    %squeezing parameter
N   = length(r);                                 %squeezed channels
M   = p.modes;                                   %total modes
t   = p.tr';                                     %transmission correction
re  = p.re';                                     %recycling amplitude
tr  = sqrt(1-re.^2);                             %recycled transmission
n   = sinh(r).^2;                                %photon number
m   = (1-p.thermal').*cosh(r).*sinh(r);          %initialize coherence
n(N+1:M) = 0;                                    %initialise all channels
m(N+1:M) = 0;                                    %initialise all channels
x   =  sqrt((n+m+sm)/2).*randn(M,ens);           %x-quadrature noise
y   =  sqrt((n-m+sm)/2).*randn(M,ens);           %y-quadrature noise
alp = tr.*(x+1i*y)+re.*a(1:M,:);                 %input amplitude noise

%ADD OUTPUT LOSSES, TRANSFORMATIONS, EXTRA QUANTUM NOISE

T   = p.matrix.*t;                               %corrected transmission

if isfield(p,'permute') == 1                     %check for  permutation
    T = T(p.permute,:);                          %apply random permutation 
end                                              %end permuation check

if p.tmss == 1                                   %Check two-mode squeezed state
    Hbs = 1/sqrt(2)*[1,-1;1,1];                  %Create beamsplitter matrix
    Hbsc = repmat({Hbs},1,N/2);                  %Multiple beamsplitters in cell 
    B = blkdiag(Hbsc{:});                        %Matrix direct sum
    sz = size(B);
    if sz(1) == N && sz(2) == N
        B(sz(1) + 1:M,sz(2) + 1:M) = 0.0;
    end 
    T = T*B;                                     %Multiply transmission and BS
end 

if  sm > 0                                       %if not normal ordered
  a = T*alp;                                     %matrix transformation     
  E   = eye(M)-T'*T;                             %difference from unitarity 
  if norm(E) > 1.e-12                            %if T matrix is nonunitary 
    [U,T] = schur(E);                            %Schur decomposition 
     E1   = U*sqrtm(T)*U';                       %square root of diagonal
     w    = (randn(M,ens)+1i*randn(M,ens));      %quantum noise
     w    = w*sqrt(sm/2);                        %quantum noise normalized
     a    = a+E1*w;                              %output alpha amplitudes
  end                                            %end if nonunitary
else                                             %end if not normal
  bet  = tr.*(x-1i*y)+re.*a(1+M:2*M,:);          %input conjugate noise
  ngen = alp.*bet;                               %generated input number
  alp  = T*alp;                                  %matrix transformation     
  bet  = conj(T)*bet;                            %conjugate transformation 
  a    = [alp;bet;ngen];                         %pack stochastic data
end                                              %end function