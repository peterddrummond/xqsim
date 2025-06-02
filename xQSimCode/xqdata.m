function av = xqdata(a,n,p)  
%   av = XQDATA(a,n,p) stores data averages or probabilities
%   Input is the 'a' cell array of stochastic variables, stored in parallel
%   with dimensions of: [internal indices, ensemble index].
%   Input index n gives the number of the observe function.
%   Input struct p gives parameters and observe function handles.
%   If p.bin{n} is not empty, stores probability density in bins.
%   Returned array 'av' is the average observable array, at current time.
%   The  bin index is used to index over probabilities, or else omitted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by xqensemble
%   Calls     observe{n}
%   xQSIM functions licensed by Peter D. Drummond, (2025) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CHECK THE n-th OBSERVE FUNCTION AND MEASURE VARIABLES IF NOT EMPTY

if isempty (p.observe{n})                        %%check if observe defined
    av = [];                                     %%no observe defined
    return;                                      %%return if no function  
end                                              %%end function check
p.nobserve = n;                                  %%Store function index
obs = p.observe{n}(a,p);                         %%Call observe function
sz = size(obs);

%%%%%%   GET PROBABILITIES IF BINS AVAILABLE
                                             
if isempty(p.binranges{n})                       %%If probabilities needed
    av = mean(obs,p.length{p.seq}{n});           %%Take average of ensemble
else
    binranges = p.binranges{n};                  %%Store the n-th bins
    or  = real(obs);                             %%Get real part of data
    k   = ones(1,sz(2));                         %%Initialize total indices
    in  = ones(1,sz(2));                         %%Initialize inrange 
    tot = 1;                                     %%Initialize total size
    m   = 0;
    for m1 = 1:min(length(binranges),sz(1))      %%Loop over non-empty bins
      pl   = length(binranges{m1})-1;
      if pl>0
        m = m+1;
        b1 = binranges{m1}(1);                   %%Get start of m-th bin
        db = binranges{m1}(2) - binranges{m1}(1);%%Get delta of m-th bin
        om = floor((or(m1,:)- b1)/db);           %%Calculate m-th index  
        inrange = (om>=0).*(om<pl);              %%Is m-th index in range?
        om = om.*inrange;                        %%Make m-th index in range
        k  = k + om*tot;                         %%Compute total indices
        tot = tot*pl;                            %%Get cumulative bin size
        in = inrange.*in/db;                     %%All indices in range?
      end
    end                                          %%End loop over variables
    in=in/p.ensembles(1);                        %%Normalise probabilities
    indsize(1) = sz(2)/p.ensembles(1);           %%Index dimension 1
    av  = zeros(indsize(1),tot);                 %%Reshape av for binning
    indsize(2) = p.ensembles(1);                 %%Index dimension 2
    k = reshape(k,indsize);                      %%Reshape index array
    in = reshape(in,indsize);                    %%Reshape range switch 
    for s = 1:indsize(1)                         %%Loop over first index
      for e = 1:indsize(2)                       %%Loop over ensemble index
        av(s,k(s,e)) = av(s,k(s,e)) + in(s,e);   %%Increment probability
      end                                        %%End loop over ensemble
    end                                          %%End loop over 1st index
    if isfield(p,'limits')&&length(p.limits)>=n&&~isempty(p.limits{n}{2})
      av = max(av, p.limits{n}{2}(1));           %%Truncate to limits
    end
end                                              %%End if probabilities
end                                              %%End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XDATA FUNCTION