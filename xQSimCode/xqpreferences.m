function in = xqpreferences(in)
% in = xqpreferences(in) sets structure preferences for xqsim
% Input and output are arrays of parameter structures 
% Licensed by Peter D. Drummond, (2025) - see License.txt 

% INITIALIZE STANDARD INPUTS

if ~iscell(in)                                        %if input not a cell
  in = {in};                                          %make input a cell
end                                                   %end if not cell
ls = length(in);                                      %get number of cells
for seq = 1:ls                                        %loop over sequence
p = in{seq};                                          %get parameters
xqvalid(seq,p);                                       %validate labels
p.name =     qprefer(p,'name',1,'');                  %Default name
p.dimensions = 0;                                     %dimensions=0 for graphs
p.seq  =     qprefer(p,'seq',1,seq);                  %Default sequence 
p.modes =    qprefer(p,'modes',1,1);                  %Default matrix size
p.minvar =   qprefer(p,'minvar',1,10^-100);           %set minimum variance
p.thermal =  qprefer(p,'thermal',p.modes,0);          %Default decoherence
p.tr    =    qprefer(p,'tr',1,1);                     %Default transmission
p.tmss =     qprefer(p,'tmss',1,0);                   %Default 2-mode squeezing
p.re    =    qprefer(p,'re',1,0);                     %Default recycling
p.observe =  qprefer(p,'observe',1,{@(a,p) a});       %Get default observe
p.pnames  =  qprefer(p,'pnames',1,{'+P ','W ','Q '}); %phase-space names
p.phase  =   qprefer(p,'phase',1,1);                  %phase-space type
p.pname  =   qprefer(p,'pname',1,p.pnames{p.phase});  %define current name
p.gen     =  qprefer(p,'gen',1,@xqgen);               %Get default gen
p.observe =  qprefer(p,'observe',1,{@(a,p) a(1,:)});  %Get default observe
p.averages = qprefer(p,'averages',1,1:length(p.observe));%Get default graphs
graphs     = max(p.averages);                         % graphs needed
p.matrix =   qprefer(p,'matrix',1,@Identity);         %Set matrix function
p.projmp  =  qprefer(p,'projmp',1,2);                 %Matrixp projections
if isnumeric(p.matrix)
  p.tname =  qprefer(p,'tname',1,'numerical');        %Default matrix name
else
  p.tname =  qprefer(p,'tname',1,p.matrix(0));        %Default matrix name 
end
p.verbose =  qprefer(p,'verbose',1,1);                %Default verbose mode
p.ensembles= qprefer(p,'ensembles',1,[1,1,1]);        %Default ensembles
p.sqz =      qprefer(p,'sqz',p.modes,0);              %Default squeezing
p.counts =   qprefer(p,'counts',1,0);                 %Default counts
p.mincount = qprefer(p,'mincount',1,10);              %set minimum count
p.binranges= qprefer(p,'binranges',1,{[]});           %set probability bins
p.compare =  qprefer(p,'compare',1,{[]});             %Default compare 
p.cyc     =  qprefer(p,'cyc',1,1);                    %Default cycle 
p.diffplot=  qprefer(p,'diffplot',1,{[]});            %Difference plot type
p.cutoff  =  qprefer(p,'cutoff',1,-1.e100);           %Default lower cutoff
p.olabels =  qprefer(p,'olabels',graphs,{' '});       %Default label;                    
p.errors  =  3;                                       %define error type
p.alpha   =  qprefer(p,'alpha',p.modes,0);            %Default coherent alpha;

% INITIALIZE NON-STANDARD INPUTS

if seq == 1                                           %if initial sequence
    rep = p.ensembles(2)*p.ensembles(3);              %get first repeats
end                                                   %end if initial
p.rep = rep;                                          %make repeats uniform
while  length(p.observe) < max(p.averages)            %if graphs too large
    p.averages =  p.averages(1:end-1);                %reduce graphs
end
graphs = max(p.averages);                             %end if graphs
p.cutoffs =  qprefer(p,'cutoffs',graphs,{p.cutoff});  %Default lower cutoff
p.binranges{graphs+1} = [];                           %extend binranges
p.initialobs{graphs+1} =[];                           %set up initial obs
p.compare{graphs+1}=[];                               %initialize compare
p.olabels{graphs+1}=[];                               %initialize labels
p.part{graphs+1}=[];                                  %initial partitions
p.xk{graphs+1}=[];
maxcount = p.modes; %initial axis labels
for k = p.averages                                    %loop over data types
    if ~isempty(p.initialobs{k})                      %if initialization?
        p = p.initialobs{k}(p);                       %initialize data type
    end                                               %end if requires
    if  isempty(p.olabels{k})                         %has no label?
        p.olabels{k} = sprintf('Graph %d',k);         %fill empty label
    end                                               %end if no label
    if  isempty(p.part{k})                            %has no partition?
        p.part{k} = {1:p.modes};                      %fill partition
    end                                               %end if no partition
    if ~iscell(p.part{k})
        vec = p.part{k};
        nd = length(vec);
        p.part{k}=cell(1,nd);
        m = 0;
        for j=1:nd
            m2 = m+vec(j);                            %add vector component
            p.part{k}{j}=(m+1):m2;
            m = m2;
        end                                           %end dimension loop
    end                                               %end if not cell
    if ~isempty(p.xk{k})                              %store max count
        nd = length(p.xk{k});
        for n = 1:nd
            maxcount = max(maxcount,p.xk{k}{n}(end));
        end
    end 
end                                                   %end loop on averages
%p.min     =  qprefer(p,'min',graphs,0);               %Default min count; 
%p.max     =  qprefer(p,'max',graphs,maxcount);        %Default max count;         
p.lfact   =  zeros(maxcount+1,1);                        %lfact(n)=log((n-1)!)
%maxcount  =  max(p.max,maxcount);
for n = 2:maxcount                                    %store factorials
    p.lfact(n+1,1) = log(n)+p.lfact(n,1);
end
if ~isnumeric(p.sqz)                                  %squeeze function
    p.sqz = p.sqz(p);                                 %get squeezing vector
end                                                   %end if from function
if ~isnumeric(p.matrix)                               %matrix from function
    p.matrix = p.matrix(p);                           %get transmission
end                                                   %end if from function
while length(p.ensembles) < 3                         %while length small
    p.ensembles(1+length(p.ensembles)) = 1;           %increment ensembles
end                                                   %end if length small
N = length(p.sqz);
if  ~isequal(p.thermal,0) 
  if length(p.thermal) < N    %check decoherence
    fprintf('Sequence %d: thermal extended\n',seq);   %vector size mismatch
    p.thermal(1,end:N) = p.thermal(1,end);            %duplicate end entry
  end                                                 %end check length
  p.thermal(1,end+1:p.modes) = 0;
  p.thermal(:,p.modes+1:end) = [];                    %reset size
end
p.sqz(N+1:p.modes) = 0;                               %reset size
p.sqz(p.modes+1:end) = [];                            %reset size
if  length(p.tr) < p.modes                            %check transmission
  fprintf('Sequence %d: tr extended\n',seq);          %vector size mismatch
  p.tr(1,end:p.modes) = p.tr(1,end);                  %duplicate end entry
end  

% INITIALIZE OUTPUT DATA AND PROBABILITY SIZES 

a = ones(3*p.modes,2);                                %initialize a
for n = p.averages                                    %loop over data types
    p.nobserve = n;                                   %store graph number
    D1 = p.observe{n}(a,p);                           %get vector average
    sz = size(D1);                                    %get data size
    if sz(end) > 2                                    %if already averaged
        sz(end+1)=1;
    else
        sz(end) = 1;
    end
    nd=length(sz)-1;
    lines = sz(1);
    p.length{seq}{n} = length(sz);                    % get data length
    p.d.bins{n} = [];                                 % Initial dimension
    binranges  = p.binranges{n};                      % n-th bin ranges
    p.bins{n} = 0;                                    % n-th bin numbers
    if ~isempty(binranges)                            % If n-th bins set         
        p.bins{n} = min(length(binranges),lines);
        dbin = zeros(1,p.bins{n});
        abin = 1;                                     % Initial area of bin
        m = 0;
        for m1 = 1:p.bins{n}                          % Loop on axes 
          d = length(binranges{m1})-1;
          if d > 0
             m = m+1;
             dbin(m) = d;                        % Dimension of n-th bin
             p.oc{n}{m} = (binranges{m1}(2:d+1)+binranges{m1}(1:d))/2;
             p.do{n}{m} = (binranges{m1}(2:d+1)-binranges{m1}(1:d));
             abin = abin*p.do{n}{m}(1);          % Area of n-th bin          
             p.xk{n}{nd+m}  = p.oc{n}{m};        % Make a new axis vector
          end
        end                                      % End loop over length
        p.bins{n} = m;                           % n-th bin numbers
        if m > 0
          if p.verbose > 0                       % If verbose printing
            fprintf ('\n%d#%d is a %d D probability',seq,n,p.bins{n});
          end                                    % end if verbose printing
          p.d.bins{n} = dbin;                    % Dimensions of bin
          p.a.bins{n} = abin;                    % Area of bin
          sz = [sz(2:end-1),dbin,sz(end)];       % Dimension of data+bins
        end
    end                                          % End if probability bins
    p.els1{seq}{n} = [1,prod(sz),1];             %get  data element
    p.els{seq}{n} =  [p.cyc,prod(sz),3];         %get total  elements
    sz(end) = 3;                                 %reset data size
    if p.cyc>1                                   %check if extra cycles
        p.sz{seq}{n} = [1,p.cyc,sz];             %get data size + cycles
    else                                         %else no extra cycles
        p.sz{seq}{n} = [1,sz];                   %get data size 
    end                                          %end check if cycles
end                                              %end loop over types

% PRINT INPUTS IF IN VERBOSE MODE
 
if p.verbose >= 2                                %if verbose  flag
    if seq == 1                                  %if initial sequence
      fprintf('\n xqpreferences v2.0 \n');       %print version name
      fprintf('\n verbose mode = %d\n',p.verbose);%verbose mode
    end                                          %end if initial 
    fprintf('\n Sequence = %d\n',seq);           %display sequence
    fprintf('\n Squeezing  \n\n') ;              %squeeze header                                
    disp(p.sqz);                                 %display squeezing  
    fprintf('\n Transmission matrix  \n\n') ;    %transmission header                                
    disp(p.matrix);                              %display matrix
    fprintf('\n Parameters \n\n') ;              %parameter head                                
    disp(p);                                     %display parameters    
end                                              %end if verbose flag    
in{seq} = p;                                     %return parameters
end                                              %end sequence loop

end                                              %end function

function pn = qprefer(struct,label,len,default)
% pn = QPREFER(struct,label,len,default) 
%   if no parameter, sets the default
%   Input: structure: 'struct', label: 'label', length: len
%   default values `default'  
%   Output: preferred value for struct.label if none is present.
%   Called by: xqpreferences
%   xqsim licensed by Peter D. Drummond, (2025) - see License 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET DEFAULTS IF NO INPUT

if ~isfield(struct,label)                        %%no input data present
    pn = default;                                %%set to default value
else
    pn = struct.(label);                         %%set to input value
    if iscell(default) 
        if ~iscell(pn)
            pn = {pn};
        end
    end
end                                              %%end if no input data
if iscell(pn) && length(pn) < len
    pn{len+1} = [];  
    for ind = 1:len
        if isempty(pn{ind})
          pn{ind} = default{1};
        end
    end
end
end