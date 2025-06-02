function [es,data,output] = xqsim(varargin)
% [error,data,output] = XQSIM(input) simulates quantum correlations in 
% photonic networks using quantum phase-space representation phases.
% Inputs come from a cell array of structures: input = {p_1,..p_k}.
% Each successive input is from a structure p of input parameters.
% Returns an error, simulation data, and output parameters with defaults
% qcpsim uses +P, W or Q  phase-space algorithms, set by phase  = 1,2,3
% If matrix  = @Unitary, uses a randomized set of unitary network matrices
% If matrix  = @Identity, uses an identity or pass-through network matrix
% Other matrix functions are allowed, as well as direct numerical inputs
% Gaussian random numbers are used to generate quantum phase-space samples
% Thermalized squeezed inputs are input to n channels of an m*m matrix
% Inputs are specified with a squeezing parameter r, decoherence thermal
% A transmission amplitude t is used to rescale the unitary
% A recycling amplitude re is used if the fields are recycled for sequences
% Calculates correlations, variances specified in s.observe functions
% For +P distributions only, it can calculate mode click probabilities. 
% Analytic comparisons or experimental data, are input in compare functions
% All output data may include error bars, accessed by the last data index
% verboseing set by verbose = 0 (terse),  = 1 (normal), = 2 (verbose)
% Version 2.0 - uses a multiply partitioned output count - see xQSim manual
% Data output compatible with the xGRAPH multidimensional graphics program
%   Licensed by Peter D. Drummond, (2021) - see License.txt 

% INITIALIZE PARAMETERS AND DATA

tic();                                           %start timer
if iscell(varargin{1})
    varargin = varargin{1};
end
fprintf('\n xqsim, v3.0 \n');                    %verbose version
output = xqpreferences(varargin);                %include default inputs
sequence = length(output);                       %get sequence length
p = output{1};                                   %get first structure
pens = p.ensembles(3); ns = p.ensembles(1);      %store ensembles
fprintf('\n%d repeats of %d samples\n',p.rep,ns);%verbose ensembles

%SIMULATE WITH PARALLEL OR SERIAL ENSEMBLE FOR DATA AVERAGES

graphs = max(p.averages);
data(1:sequence) =  {cell(graphs,1)};            %initialize graph cell
if pens <= 1                                     %if no parallel ensembles
    data = xqensemble(1,output);                 %compute ensemble data
else                                             %else parallel ensembles
    parfor r = 1:pens                            %do parallel loop
        data1 = xqensemble(r,output);            %compute ensemble data
        data = qaddcell(data,data1);             %average ensemble data
    end                                          %end parallel loop
end                                              %end if no parallel

% LOOP OVER DATA TYPE TO GET STANDARD DEVIATIONS

es = 0; ed = 0; chi = 0.; ecmp = 0; chn = 0;     %initialize errors                            
r1 = p.rep-1;                                    %reduce for errors in mean

for s = 1:sequence                               %loop over sequence
  p = output{s};                                 %parameters for sequence
  fprintf('\nDataset %d, %s phase-space \n',s, p.pname); % sequence number
  fprintf('%d x %d %s matrix\n',p.modes,p.modes,p.tname);%matrix used
  fprintf('Squeezed input = %d\n',length(p.sqz));%input size
  for n = p.averages                             %loop over data type
    esk = 0; edk = 0; ecmpk = 0;                 %initialize k errors
    fprintf('\nDataset %d, Graph %d, %s\n',s,n,p.olabels{n}); % k,s
    sz = size(data{s}{n});                       %get size of the data
    els = prod(sz(1:end-1));                     %length of the data type
    data{s}{n} = reshape(data{s}{n},[els,3]);    %reshape data for summing
    av = data{s}{n}(:,1);                        %get averages
    if r1>0                                      %if standard deviations
      sd  = sqrt((data{s}{n}(:,3) - av.^2)/r1);  %get standard deviations
      esk = max(esk,max(real(sd)));              %store kmax. standard dev.
      es = max(es,esk);                          %store max. standard dev.
      data{s}{n}(:,3)   = sd;                    %store the standard dev.
      fprintf('\nMax samp. error = %d \n',esk);  % k sampling error
    else                                         %no standard deviations
      sd = 0;                                    %make standard dev. = 0
      data{s}{n}(:,3)   = 0;                     %store standard dev. = 0
      fprintf('\nNo sampling error available\n');%verbose warning
    end                                          %end if standard dev.
    bin=p.binranges{n};
    if ~isempty(bin)   
      output{s}.xk{n} = output{s}.xk{n}(2:end);  %shift xk{n}: probability.
      db = 1;
      for m1 = 1:length(bin)
        db = db*(bin{m1}(2) - bin{m1}(1));    %%Get delta of m-th bin
      end
      prob = sum(av)*db;
      fprintf('\nOutput = %d,total probability = %d, \n',n,prob);
    end

% COMPARISONS AND CHI SQUARE TEST
    
    if ~isempty(p.compare{n})                    %check if comparisons at k
        p.noutput = n;                           %store graph number
        cp = reshape(p.compare{n}(p),[els,1]);   %reshape comparisons
        cp = cat(2,cp,zeros(els,2));             %add error arrays
        del = abs(cp(:,1)- av);                  %get differences
        edk = max(edk,max(del));                 %store max. difference
        ecmpk = max(ecmpk,max(cp(:,2)));         %store max. compare error
        ed = max(ed,edk);                        %store max. difference
        ecmp = max(ecmp,ecmpk);                  %store max. compare error
        nz = cp(:,1) > p.cutoffs{n};             %cutoff switch
        if p.counts > 0                          %if counts exist
          nz = nz & (p.counts*av>=p.mincount);   %update cutoff switch
          cp(:,3) = sqrt(av/p.counts);           %get estimated sd
        end                                      %end if counts exist
        var = sd.^2 + cp(:,3).^2 + p.minvar;     %combined variances
        nzk = sum(nz); chn = chn+nzk;            %k significant points
        chk = sum(nz.*(del.^2./var));            %chi-squared data for k
        chi = chi+chk;                           %all chi-squared totals        
        fprintf('Max diff. error = %d \n',edk);  % difference error
        fprintf('Max comp. error = %d \n',ecmpk);% comparison error
        if r1>0 
        fprintf('Chi-square points = %d \n',nzk);% comparison points
        fprintf('Chi-square error = %d\n',chk);  % k chi-square error
          if nzk >= 10                         
            whchk = (chk/nzk)^(1/3);             % wilson-hilferty(WH) 
            whvar = 2/(9*nzk);                   % WH variance
            Zst = (whchk-(1-whvar))./sqrt(whvar);%compute WH Z-stat
            fprintf('WH Z-statistic = %d\n',Zst);%print WH Z-statistic
          end
        end
        data{s}{n}  = cat(2,data{s}{n},cp);      %combine data for xgraph;
        sz(end)  = 6;                            %update data size
    end                                          %end check comparison
    data{s}{n}  = reshape(data{s}{n},sz);        %reshape data for xgraph;
  end                                            %end loop over data
end                                              %end loop over sequence
fprintf('\nSummary of all %d xqsim datasets\n',sequence);%total error
fprintf('Max sampling error = %d \n',es);        %verbose total error
if chi > 0 ||  ed > 0  || ecmp > 0
  chi = chi/(chn+1.e-100);                       %get average CS error
  fprintf('Max difference error = %d \n',ed);    %difference error
  fprintf('Max comparison error = %d \n',ecmp);  %comparison error
  if r1>0 
    fprintf('Mean chi-square error = %d \n',chi);%average CS error
     es = chi+es*(chi==0);                       %return chi if nonzero
  end
end
fprintf('Simulation time = %f \n\n\n',toc());    %time taken
end                                              %end simulation function

% DEFINE FUNCTION TO ADD DATA CELLS

function d = qaddcell(d,data)
%   d = QADDCELL(d,data) 
%   Adds two cell arrays together, initialize first if empty
%   Licensed by Peter D. Drummond, (2021) - see License.txt 

for seq = 1:length(data)
  for n = 1:length(data{seq})
      if isempty(d{seq}{n})
        d{seq}{n} = zeros(size(data{seq}{n}));  %initialize if empty
      end
      d{seq}{n} = d{seq}{n}+data{seq}{n};       %add to returned cells
  end
end
end
