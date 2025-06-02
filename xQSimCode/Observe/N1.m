function C = N1(a,varargin)
% C = N1(a,p); generates count probabilities in a single partition
% Input is EITHER N1(a,p) or N1(a,modes,counts,p)
% If input: modevec is a vector of modes included in the sum
%           counts is a vector of counts to compute probability
% Otherwise modevec  = p.part{p.nobserve}{1};
% Also,     counts   = p.xk{p.nobserve}{1};
% Returns probability of counts on list for all modes counted
% This is given by: P(c) = (1/c!)(np^c)exp(-np)
% First output index = c1 = n+1, second index = ensembles

switch nargin
  case 4
    modevec  = varargin{1};
    counts   = varargin{2};
    p        = varargin{3};
  case 2
    p        = varargin{1};
    n        = p.nobserve;
    modevec  = p.part{n}{1};                     
    if isempty(p.xk{n}) || isempty(p.xk{n}{1})   %if no count numbers
        counts = 0:p.modes;
    else
        counts = p.xk{n}{1};
    end
  otherwise
    error('N1 requires 2 or 4 inputs, not %d', nargin);
end
if p.phase ~= 1
    error('N1 requires +P method = 1');
end
np = sum(a(modevec,:).*a(p.modes+modevec,:),1);  % photon numbers
lp = log(np);                                    % log photon numbers
C  = exp(counts'.*lp-np-p.lfact(1+counts,1));    % count probabilities   
end                                              % end N1 function