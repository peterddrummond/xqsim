    function input = xgpreferences (input,oldinput)  
%   input  =  XGPREFERENCES(input) sets default values for  graphics inputs.
%   Input:  'input' cell array, and any previous input 'oldinput'.
%   Output: 'input' cell array with updated and default graphics values.
%   Called by: xgraph
%   Needs:     xprefer, xcprefer
%   xGRAPH functions are licensed by Peter D. Drummond, (2021) - see License
%
fprintf ('\nxGRAPH, v3.4\n');                       %%version name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  COMPARE NEW AND OLD SEQUENCES
%                          
if ~isempty(oldinput)                            %% Check  oldinput ~empty
    oldsequence = length(oldinput);              %%get sequence length
    for s= 1:oldsequence                         %%loop over old sequence 
        in = oldinput{s};                        %%get old sequence input
        fname = fieldnames(in);                  %%get old sequence  labels
        if s<= length(input)                     %%check new input length
          for j = 1:length(fname)                %%loop over old labels
            label = fname{j} ;                   %%set new label = old 
            if ~isfield(input{s},label)          %% If no label data 
              input{s}.(label) = in.(label);     %% Set new input to old
            end                                  %% End if no label data
          end                                    %% End for loop
        else                                     %%check no new input
          input{s} = oldinput{s};                %%set new input to old
        end                                      %% End s < = length
    end                                          %% End sequence loop
end                                              %% End if oldinput ~empty
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET THE STANDARD DEFAULT VALUE
%                          
%input = xpreferences(input);                     %%get any default inputs                          
sequence = length(input);                        %%get sequence length
linestyle = {'-k','--k',':k','-.k','-ok','--ok',':ok','-.ok','-+k','--+k'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE WHOLE SEQUENCE
%
dgraphs = input{1}.dgraphs; 
for s = 1:sequence                               %%loop over sequence 
  in = input{s};                                 %%current inputs
  octave = exist('OCTAVE_VERSION','builtin');    %%finds if octave version
  if ~isfield(in,'dimension')                    %% If no label data
      in.dimension = 0;
  end
  if ~isfield(in,'numberaxis')                    %% If no label data
      in.numberaxis = 0;
  end
  if ~isfield(in,'functions')                    %% If no functions data
      in.functions = 0;
  end
  maxd=10;                                       %%max number of dimensions 
  xlabels=cell(1,maxd);                          %%cell array for xlabels 
  klabels=cell(1,maxd);                          %%cell array for klabels
  bounds{1} = num2cell(zeros(1,maxd));           %%initialize bounds
  if in.numberaxis || in.dimension > 4
      for i=1:in.dimension
          xlabels{i} = sprintf('r_%d',i);
          klabels{i} = sprintf('k_%d',i);
       end 
  else
      if in.dimension
        xlabels={'t','x','y','z'};
        klabels={'\omega','k_x','k_y','k_z'};
      else
        xlabels={'m_1','m_2','m_3','m_4','m_5','m_6'};
      end
  end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET COMMON GRAPHICS DEFAULTS
%                          
  in.octave =     xprefer(in,'octave',1,octave,0,5);
  in.print  =     xprefer(in,'print',1,1,0,0);
  in.graphs =     xprefer(in,'graphs',1,dgraphs,0,0);
  in.gtransforms= xcprefer(in,'ftransforms',in.graphs,{zeros(1,maxd)});
  in.minbar =     xcprefer(in,'minbar',in.graphs,{0.01});
  in.cutoff =     xprefer(in,'cutoff',1,-1.e100,0,0);
  in.cutoffs =    xcprefer(in,'cutoffs',in.graphs,{in.cutoff});
  in.linestyle =  xcprefer(in,'linestyle',in.graphs,{linestyle});
  in.legends =    xcprefer(in,'legends',in.graphs,{''});
  in.limits   =   xcprefer(in,'limits',in.graphs,{zeros(maxd,2)});
  in.esample =    xcprefer(in,'esample',in.graphs,{1});
  in.font =       xcprefer(in,'font',in.graphs,{18});
  in.headers =    xcprefer(in,'headers',in.graphs,{''});
  in.logs    =    xcprefer(in,'logs',in.graphs,{zeros(1,maxd)});
  in.images =     xcprefer(in,'images',in.graphs,{0});  
  in.imagetype =  xcprefer(in,'imagetype',in.graphs,{1});
  in.slices =      xcprefer(in,'slices',in.graphs,{0});
  in.transverse = xcprefer(in,'transverse',in.graphs,in.slices);
  in.bounds     = xcprefer(in,'bounds',in.graphs,bounds); 
  in.pdimension = xcprefer(in,'pdimension',in.graphs,{3});
  in.diffplot = xcprefer(in,'diffplot',in.graphs,{0});
  in.compares =   xcprefer(in,'compares',in.graphs,{0});
  in.parametric = xcprefer(in,'parametric',in.graphs,{zeros(1,2)});
  in.compare = xcprefer(in,'compare',in.graphs,{''});
  in.olabels = xcprefer(in,'olabels',in.graphs,{' '});
  in.errors =   xprefer(in,'errors',1,1,0,0); %%Number of error fields
  if ~isfield(in,'gfunction') && in.graphs > 0     %% If no label data
      for n = 1:in.graphs
        in.gfunction{n} = @(d,~) d;
      end
  end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE GRAPH INDEX
%
  for n = 1:in.graphs                            %% Loop over graphs
    if sequence == 1 
        stn = sprintf(' #%d ',n);
    else
        stn = sprintf(' %d#%d ',s,n);
    end
    in.headers{n} = [in.name,stn,in.headers{n}]; %% Add name to header
    in.xfunctions{n}{maxd}=[];
    in.glabels{n}{maxd}=[];
    in.xk{n}{maxd} = [];
    if isfield(in,'xlabels')
        xlabels = in.xlabels;
    end
    if isfield(in,'klabels')
        klabels = in.klabels;
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE DIMENSION INDEX
%                          
    for nd = 1:maxd                              %% Loop over dimension  
      if isempty(in.glabels{n}{nd}) && nd<=length(xlabels)
           if nd<=in.dimension && in.gtransforms{n}(nd)
              in.glabels{n}{nd}=klabels{nd};
           else
              in.glabels{n}{nd}=xlabels{nd};
           end
      end                                        %% End if undefined
      if isempty(in.xfunctions{n}{nd})
            in.xfunctions{n}{nd} = @(x,~) x;     %% Return default
      end
    end                                          %% End loop over dimension
  end                                            %% End loop over graphs
  input{s} = in;                                 %% return cells for output
end                                              %% end sequence loop
end                                              %% end function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XGPREFERENCES