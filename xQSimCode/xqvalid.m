function [] = xqvalid(s,p)
 %   e = XQVALID(s,p) validates an xqsim parameter structure
 %   Here s is the sequence index and p is the parameter struct
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 validinput = {...     
  ...                   % VALID XQSIM PARAMETER LABELS
 'alpha',...
 'averages',...
 'binranges',...
 'compare',...
 'correl',... 
 'cutoff',...
 'cutoffs',...
 'cyc',...
 'deco',...  
 'ensembles',... 
 'file',...
 'gen',...
 'matrix',...  
 'mincount',...
 'modes',... 
 'name',...
 'observe',... 
 'octave',...
 'part',...
 'permute',...
 'pnames',...
 'phase',...
 're',...
 'seed',...
 'sqz',... 
 'thermal',...
 'tname',...
 'tr',...  
 'verbose',... 
 'xk',...
  ...                % VALID XGRAPH PARAMETER LABELS
 'axes',... 
 'chisqplot',... 
 'cutoffs',... 
 'diffplot',... 
 'font',... 
 'gfunction',... 
 'glabels',... 
 'graphs',... 
 'gsqplot',...
 'headers',... 
 'images',... 
 'imagetype',... 
 'legends',... 
 'limits',... 
 'linestyle',... 
 'linewidth',... 
 'logs',...
 'minbar',... 
 'mincount',...
 'olabels',...  
 'savefig',... 
 'savegraph',...
 'scale' 
}; 

pfields = fieldnames(p);
valid   = ismember(pfields,validinput);
if ~all(valid)
  OK = 1;
  for n = 1:length(pfields)
    if ~valid(n) && isstrprop(pfields{n}(1),'lower')
      fprintf('\n Invalid parameter: %s\n',pfields{n});
      OK = 0;
    end
  end
  if ~OK
    error ('\n Invalid input structure, sequence %d\n',s);
  end
end
end