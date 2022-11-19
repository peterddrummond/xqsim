function in = xqpreferences(in)
% in = xqpreferences(in) sets structure preferences for qcpsim
% Input and output are arrays of parameter structures 
% Licensed by Peter D. Drummond, (2021) - see License.txt 

% INITIALIZE STANDARD INPUTS

if ~iscell(in)                                      %if input not a cell
  in={in};                                          %make input into a cell
end                                                 %end if input not cell
ls = length(in);                                    %get number of cells
for s = 1:ls                                        %loop over sequence
p = in{s};                                          %get parameters
p.name =     qprefer(p,'name','');                  %Default name
p.seq  =     qprefer(p,'seq',s);                    %Default sequence index
p.cyc  =     qprefer(p,'cyc',1);                    %Default cycles
p.M =        qprefer(p,'M',1);                      %Default matrix size
p.N =        min(p.M,qprefer(p,'N',p.M));           %Default input size
p.minvar =   qprefer(p,'minvar',10^-100);           %set minimum variance
p.eps    =   qprefer(p,'eps',0);                    %Default decoherence
p.t    =     qprefer(p,'t',1);                      %Default transmission
p.re    =    qprefer(p,'re',0);                     %Default recycling
p.observe =  qprefer(p,'observe',{@(a,p) a});       %Get default observe
p.pnames  =  qprefer(p,'pnames',{'+P ','W ','Q '}); %phase-space names
p.method  =  qprefer(p,'method',1);                 %phase-space method
p.pname  =   qprefer(p,'pname',p.pnames{p.method}); %define current name
p.gen     =  qprefer(p,'gen',@xqgen);                %Get default gen
p.observe =  qprefer(p,'observe',{@(a,p) a(1,:)});  %Get default observe
p.graphs =   qprefer(p,'graphs',length(p.observe)); %Get default graphs
p.matrix =   qprefer(p,'matrix',@Identity);         %Set matrix function
p.Tname =    qprefer(p,'Tname',p.matrix(0));        %Default matrix name
p.print =    qprefer(p,'print',1);                  %Default print mode
p.ensembles= qprefer(p,'ensembles',[1,1,1]);        %Default ensembles
p.r =        qprefer(p,'r',0);                      %Default squeezing
p.counts =   qprefer(p,'counts',0);                 %Default counts
p.mincount = qprefer(p,'mincount',10);              %set minimum count
p.compare =  qprefer(p,'compare',{[]});             %Default compare handle
p.diffplot=  qprefer(p,'diffplot',{[]});            %Difference plot type
p.cutoff  =  qprefer(p,'cutoff',-1.e100);           %Default lower cutoff
p.errors  =  3;                                     %define error type

% INITIALIZE NON-STANDARD INPUTS

if s == 1                                           %if initial sequence
    rep = p.ensembles(2)*p.ensembles(3);            %get first repeats
end                                                 %end if initial
fprintf('\nDataset %d: %s\n',s,p.name);             %print dataset name
p.rep = rep;                                        %make repeats uniform
if  length(p.observe) < p.graphs                    %if graphs too large
    p.graphs =   length(p.observe);                 %reduce graphs
end                                                 %end if graphs
p.initialobs{p.graphs+1} =[];                       %set up initial obs
p.compare{p.graphs+1}=[];                           %initialize comparisons
p.olabels{p.graphs+1}=[];                           %initialize labels
for k = 1:p.graphs                                  %loop over data types
    if ~isempty(p.initialobs{k})                    %requires initialization?
        p = p.initialobs{k}(p);                     %initialize  data type
    end                                             %end if requires
    if  isempty(p.olabels{k})                       %has no label?
        p.olabels{k} = sprintf('Graph %d',k);       %fill empty label
    end                                             %end if no label
end                                                 %end loop on data types
if ~isnumeric(p.r)                                  %squeeze from function
    p.r = p.r(p);                                   %get squeezing vector
end                                                 %end if from function
if ~isnumeric(p.matrix)                             %matrix from function
    p.matrix = p.matrix(p);                         %get transmission
end                                                 %end if from function
while length(p.ensembles) < 3                       %while length too small
    p.ensembles(1+length(p.ensembles)) = 1;         %increment ensembles
end                                                 %end while length small
if  length(p.r) < p.N                               %check squeezing length
  fprintf('Dataset %d warning, r changed \n',s);    %input size mismatch
  p.r(1,end:p.N) = p.r(1,end);                      %duplicate end entry
end                                                 %end check length
p.r(p.N+1:p.M) = 0;                                 %reset size
p.r(p.M+1:end) = [];                                %reset size
if  length(p.eps) < p.N                             %check decoherence
  fprintf('Dataset %d warning, eps changed\n',s);   %vector size mismatch
  p.eps(1,end:p.N) = p.eps(1,end);                  %duplicate end entry
end                                                 %end check length
p.eps(p.N+1:p.M) = 0;                               %reset size
p.eps(p.M+1:end) = [];                              %reset size

% INITIALIZE DATA SIZES

a = zeros(2*p.M,p.ensembles(1));                    %initialize a
for k = 1:p.graphs                                  %loop over data types
    p.k = k;                                        %store graph number
    D1 = p.observe{k}(a,p);                         %get vector average
    sz = size(D1);                                  %get data size          
    lastk = length(sz);                             %get data length
    if lastk > 2 || sz(lastk) > 1                   %check last data index
          lastk=lastk+1;                            %increment data size
    end                                             %end check last index
    p.els1{s}{k} = [1,prod(sz),1];                  %get  data element
    p.els{s}{k} =  [p.cyc,prod(sz),3];              %get total  elements
    sz(lastk) = 3;                                  %reset data size
    if p.cyc>1                                      %check if extra cycles
         p.sz{s}{k} = [1,p.cyc,sz];                 %get data size + cycles
    else                                            %else no extra cycles
         p.sz{s}{k} = [1,sz];                       %get data size 
    end                                             %end check if cycles
end                                                 %end loop over types

% PRINT INPUTS IF IN VERBOSE MODE
 
if p.print >=2                                      %if verbose print flag
    if s == 1                                       %if initial sequence
      fprintf('\n xqpreferences v1.0 \n');           %print version name
      fprintf('\n Print mode = %d\n',p.print);      %print mode
    end                                             %end if initial 
    fprintf('\n Sequence = %d\n',s);                %display sequence
    fprintf('\n Squeezing r  \n\n') ;               %squeeze header                                
    disp(p.r);                                      %display squeezing  
    fprintf('\n Transmission matrix  \n\n') ;       %transmission header                                
    disp(p.matrix);                                 %display matrix
    fprintf('\n Parameters \n\n') ;                 %parameter head                                
    disp(p);                                        %display parameters    
end                                                 %end if print flag    
in{s} = p;                                          %return parameters
end                                                 %end sequence loop
end                                                 %end function

function pn = qprefer(struct,label,default)
% pn = QPREFER(struct,label,default) if no parameter, sets the default
%   Input: structure 'struct', label 'label', default values `default'  
%   Output: preferred value for struct.label if none is present.
%   Called by: xqpreferences
%   QSIM functions are licensed by Peter D. Drummond, (2021) - see License 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET DEFAULTS IF NO INPUT
%
   if ~isfield(struct,label)                        %%no input data present
      pn = default;                                 %%set default value
   else
     pn = struct.(label);                           %%set default value
   end                                              %%end if no input data
end