function C = Mn1(a,p)
% C = MN1(a,p); generates count probabilities in a single partition
% Uses coherent projection matrix-P method for computing observables
% Requires a pure squeezed-state input for validity
% Needs:    modev    = p.part{p.nobserve}{1};
% Also,     counts   = p.xk{p.nobserve}{1};
% Returns probability of counts on list for all modes counted
% Output is: P(m) = g(n,m)(1/m!)(np^m)exp(-np)
% where    g(n,m)= (1-(2mod(m,2)-1)exp(2(n-np))/(1+exp(-2np)).
% Here     np is the partition photon number, n the input number
% First output index = m1 = m+1, second index = ensembles
  
if ~isequal(p.phase,1)                           %check phase-space is +P
    error('Mn1 requires p.phase = 1, not %d',p.phase);
end
if ~(isequal(p.thermal,0)||isequal(p.thermal, zeros(1,p.modes)))                               
    error('Mn1 requires p.thermal = 0, not %d',p.thermal(1)');
end
if ~isequal(p.alpha,0)                           %check no coheemt input
    error('Mn1 requires p.alpha = 0, not %d',p.alpha(1)');
end
obs = p.nobserve;
xkcell = p.xk{obs};                              %count numbers for graph
if isempty(xkcell) || isempty(xkcell{1})         %if no count numbers input
    counts = 0:p.modes;
else
    counts = xkcell{1};
end
modev  = p.part{obs}{1};
n   = sum(a(1+2*p.modes:3*p.modes,:),1);         %input boson numbers
np  = sum(a(modev,:).*a(p.modes+modev,:),1);     % photon numbers
lp  = log(np);                                   % log photon numbers
par = 1-2*mod(counts',2);
C   = exp(counts'.*lp-np-p.lfact(1+counts,1));   % count probabilities
C   = C.*((1+par.*exp(2*(np-n)))./(1+exp(-2*n)));
end                                              % end MN1 function