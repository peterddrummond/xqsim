function   xgsq(n,d,g) 
%   XGSQ(n,d,g) plots g-squared max likelihood plot data 
%   for a single xGRAPH plot function, if probability data exists.
%   Input: graph number 'n', data 'd', graph parameters 'g'.
%   Output: graphs of g-squared data and k values if gsqplot = 1
%   or Wilson-Hilferty z data if gsqplot = 2
%   Sums over all the data points in a sampled distribution
%   Results are graphed against any space-time dimensions
%   Called by: xgraph
%   Needs:     xmultigraph
%   Licensed by Peter D. Drummond, (2021) - see License.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SUM OVER CHI-SQUARE FOR PLOTS
% 
dimension = 0;                                   %%initialize the dimension
if isfield(g,'dimensions')                       %%check if dimension field
    dimension = g.dimensions;                    %%store the dimension
end                                              %%end check if field
if isequal(dimension,0)                          %%if a plot not possible?
    return                                       %%return to calling code
end                                              %%end if plot
sz       = size(d);                              %% data size
axsize   = length(sz) - 2;                       %%total axes available
nsumd    = (axsize-dimension+2):(axsize+1);      %%indices for summing
nplotd   = 1:(axsize-dimension+1);               %%indices for plotting
sz1      = prod(sz(nplotd));                     %%total plot size
sz2      = prod(sz(nsumd));                      %%total summing size
d        = reshape(d,[sz1,sz2,6]);               %%reshape count data
count    = d(:,:,1)*g.scale{n};                  %%scale simulated counts
expc     = d(:,:,4)*g.scale{n};                  %%scale expected counts
nc       = (count>g.mincount)&(expc>g.mincount); %%nonzero counts
count    = count.*nc;
expc     = expc.*nc;
expc     = expc.*(sum(count,2)./sum(expc,2));
valid    = sum(nc(:,:,1),2);                     %%get number of points
ratio    = count./(expc+1.e-100);
ratio    = (count == 0.0)+ratio;
gsqu     = 2*nc.*count.*log(ratio);              %%get chi-squared data
gsqu     = sum(gsqu,2);                          %%sum over probabilities
if g.gsqplot{n} == 2                             %%If second plot type?
    gsqu = ((gsqu./valid).^(1/3))-(1-2./(9*valid));
    gsqu = gsqu.*sqrt(valid.*4.5);               %%normalise z data
end                                              %%End if plot type
gsqu(:,:,2) = 0;                                 %%set error bars to zero
gsqu(:,:,3) = 0;                                 %%sampling errors to zero
gsqu = reshape(gsqu,[sz(nplotd),3]);             %%reshape chi-square data
if g.gsqplot{n} == 1                             %%First plot type?
    g.olabels{n}= 'G^2, k';                      %%make label
    valid(:,:,2) = 0;                            %%set error bars to zero
    valid(:,:,3) = 0;                            %%sampling errors to zero
    valid = reshape(valid,[sz(nplotd),3]);       %%reshape counts data
    gsqu = [gsqu;valid];                         %%include count data
else                                             %%Second plot type?
    g.olabels{n}='z_{G^2}';                      %%make label
end                                              %%End if plot type
xmultigraph(n,gsqu,0,g);                         %%make g-square plots
end                                              %%end function call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    END XGSQPLOT FUNCTION