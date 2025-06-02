function [] = xpr(m,p,varargin)
%   inlabel  =  Xverbose(sw,p,varargin) verboses xSPDE messages
%   verboses output if m <= p.verbose.
%   Minimal if p.verbose = -1       verboses start-up and faults
%   Brief if   p.verbose =  0:     +verboses total integration errors
%   Normal if  p.verbose =  1:     +verboses function errors and progress
%   Verbose if p.verbose =  2:     +verboses internal parameters as well
%   
%   xSPDE functions licensed by Peter D. Drummond, (2021) - see License  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  if m <= p.verbose
    fprintf (varargin{:});
  end
end      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XPR
%
