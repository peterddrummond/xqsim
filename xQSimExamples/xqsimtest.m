function [] = xqsimtest()
%   XQSIMTEST is a  script to test the operation of the xQSIM toolbox 
%   There are currently 6 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for this test.
%   If not available, replace parallel ensembles by serial ensembles
%   Runtime: 20-100s for Matlab R2024, depending on cores available
%            <1200s for Octave v6.2
%   To see the generated graphs, add '%' before the  'close all' commands
%   If needed, you can also include: 
%cd ~               %Use home directory for stored data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=tic;
e{1} = Multi500;         n{1} = 'Multi500';
e{2} = GBS8PNRMP;        n{2} = 'GBS8PNRMP';
e{3} = GBS20P;           n{3} = 'GBS20P';
e{4} = GBS20PNRP;        n{4} = 'GBS20PNRP';
e{5} = GBS20;            n{5} = 'GBS20';
e{6} = xquad1Dmarginal;  n{6} = 'xquad1Dmarginal';
close all;               %%Deletes all figures if not wanted
Et = toc(t1);      
                         %%Print summarized error scores and time
                   

fprintf('\n\nxQSIM3.0 Batchtest  errors, typical in brackets:\n\n');
fprintf(['Result: (Typical) Chi-2/k  Name\n\n']);
e1 = [2.140498e-01,6.638778e-01,9.105589e-01,4.613130e-01,1.378932e+00,...
      9.014873e-01];
Exp = 0;
et = 0;
for j = 1:length(e)
  if e{j}>e1(j)/10 && e{j}<e1(j)*10 
    fprintf('# %2d OK (%.6e) %.6e  (%s)\n',...
        j,e1(j),e{j},n{j});
  else
    fprintf('# %2d ?? (%.6e) %.6e (%s)\n',...
        j,e1(j),e{j},n{j});
  end
  et = et+e{j};
  Exp = Exp+e1(j);
end
if (et>Exp/2 && et<2*Exp) 
 fprintf('\nTotal error within 50%% of expected error\n');
else
  if (et>Exp/10 && et<10*Exp) 
    fprintf('\nTotal error within order of magnitude of expected value\n');
    fprintf('If Octave, errors may differ from Matlab error estimates \n');
  else
    fprintf('\nBatchtest failed, check installation\n');
  end
end
fprintf('\nTotal elapsed time  = %.3g (<~120) seconds\n',Et);
end