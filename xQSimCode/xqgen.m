function a = xqgen(a,p)
% a  = xqgen(a,p) generates phase-space data for boson sampling
% Input is a complex matrix, a = alpi;beti and a parameter struct p
% First the input amplitude is scaled by the recycling factor p.re
% This is combined with a local squeezed state or thermal excitation
% Squeezed state parameter is p.r, decoherence parameter is p.eps 
% p.matrix is the transmission, p.method is the phase-space method = 1,2,3 
% Transmission is modified by transmission correction factor p.t 
% Output is a = alp;bet, including decoherence and transmission factors
% First a matrix index gives the channel count
% Second a matrix index gives the index of the parallel ensemble
% *Licensed by Peter D. Drummond & Alexander S. Dellios, (2023) 
% - see License.txt, xQSim manual  

%GENERATE PHASE-SPACE SAMPLES FOR ANY ORDERING

ens = p.ensembles(1);                            %ensemble number
sm  = (p.method-1)/2;                            %ordering parameter
r   = p.r';                                      %squeezing parameter
t   = p.t';                                      %transmission correction
re  = p.re';                                     %recycling amplitude
tr  = sqrt(1-re.^2);                             %input field transmission
n   = sinh(r).^2;                                %photon number
m   = (1-p.eps').*cosh(r).*sinh(r);              %initialize coherence
n(p.N+1:p.M) = 0;                                %initialise all channels
m(p.N+1:p.M) = 0;                                %initialise all channels
x   =  sqrt((n+m+sm)/2).*randn(p.M,ens);         %x-quadrature noise
y   =  sqrt((n-m+sm)/2).*randn(p.M,ens);         %y-quadrature noise
alpi= tr.*(x+1i*y)+re.*a(1:p.M,:);               %input amplitude noise
beti= tr.*(x-1i*y)+re.*a(1+p.M:2*p.M,:);         %input conjugate noise

%ADD OUTPUT LOSSES, TRANSFORMATIONS, EXTRA QUANTUM NOISE

T   = p.matrix.*t;                               %corrected transmission

if isfield(p,'permute') == 1                     %check for  permutation
    T = T(p.permute,:);                          %apply random permutation 
end                                              %end permuation check


alp = T*alpi;                                    %matrix transformation                        
bet = conj(T)*beti;                              %conjugate transformation 
if  sm > 0                                       %if not normal ordered
E   = eye(p.M)-T'*T;                             %difference from unitarity 
  if norm(E) > 1.e-12                            %if T matrix is nonunitary 
    [U,T] = schur(E);                            %Schur decomposition 
     E1   = U*sqrtm(T)*U';                       %square root of diagonal
     w    = (randn(p.M,ens)+1i*randn(p.M,ens));  %quantum noise
     w    = w*sqrt(sm/2);                        %quantum noise normalized
     alp  = alp+E1*w;                            %output alpha amplitudes
     bet  = conj(alp);                           %output beta amplitudes
  end                                            %end if nonunitary
end                                              %end if not normal
a   = [alp;bet];                                 %pack stochastic data
end                                              %end function