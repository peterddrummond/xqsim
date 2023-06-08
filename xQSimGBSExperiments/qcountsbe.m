function [valid,prob1] = qcountsbe(p)
%[C] = QCOUNTSBE()
%   Gets the experimental count probability for a GBS experiment
%   Input data must come from big endian binary files
%   Data ends in .bin, last data point is checked for a validity flag
%   Mask is applied to remove time-stamps and specific invalid data points
%   Input data should be saved in the working directory used by Matlab.
%   Produces individual data files for each observable. 
%   Licensed by Peter D Drummond, (2020).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Removes and relabels invalid data, see data description
% V3.0
%

tic()
fprintf('ExpCounts1_0\n');                  % print current version name
maxrecords = 60000000;
%maxrecords = 600000;
fprintf('Processing up to %d records\n', maxrecords); 
%maxrecords = 600;

% 
% Mask data for big endian ordering. This section is user specified for
% each data set. 
%
mask = zeros(p.bits,1);
index  = mask;
mask(17:120) = 1;
mask(23) = 0;
mask(25) = 0;
mask(26) = 0;
mask (44) = 0;
mask(128) = 1;
for j=1:120
    index(j) = (p.M + 1)-sum(mask(1:j));
end

%
% Create observables.
%
count1 = zeros(p.M+1,1);
ln = length(p.Part);
for i = 1:ln
    sz = size(p.Part{i});
    if sz(2) == 2
        count2 = zeros(p.Part{i}+1);
    elseif sz(2) == 4
        count4 = zeros(p.Part{i}+1);
    end 
end 
click = zeros(p.M,1);
n = zeros(p.M,1);
K = 5;
Kc = zeros(K,1);

%
% Load in data from available bin files
%

d     = dir('*.bin');                       % gets files labelled *.bin
if isempty(d) 
    error('Found 0 binary files in directory\n');
end
if length(d)>1
    error('Found %d binary files\n', length(d));
else
    fprintf('Found %d binary file\n', length(d)); 
end
id = fopen(d(1).name);                     % Open bin file
bits = p.bits;
record = 0;
valid = 0;

while bits  == p.bits && record < maxrecords
    [B1,bits] = fread(id,p.bits,'ubit1','ieee-be');
    if bits == p.bits
      record =record+1;
      if B1(p.bits,1) == 0
        click(1:76)   = B1(120:-1:45);%%76
        click(77:93)  = B1(43:-1:27); %%17
        click(94)     = B1(24);       %%1
        click(95:100) = B1(22:-1:17); %%6
        n             = n+click;
        for i = 1:K
            Kc(i)     = Kc(i)+ prod(click(1:i));
        end
        c             = 1+sum(click);
        c21           = 1+sum(click(1:(0.5*p.M)));
        c22           = 1+sum(click((0.5*p.M+1):p.M));
        c41           = 1+sum(click(1:(0.25*p.M)));
        c42           = 1+sum(click((0.25*p.M+1):(0.5*p.M)));
        c43           = 1+sum(click((0.5*p.M+1):(0.75*p.M)));
        c44           = 1+sum(click((0.75*p.M+1):p.M));
        count1(c,1) = count1(c,1)+1; 
        count2(c21,c22) = count2(c21,c22)+1;
        count4(c41,c42,c43,c44) = count4(c41,c42,c43,c44)+1;
        valid = valid+1;
      end
    end
end    
fprintf(' Total records = %d\n',record);
fprintf(' Valid records = %d\n',valid); 

n    = n/valid;
Kc   = Kc/valid;
prob1 = count1/valid;
prob2 = count2/valid;
prob4 = count4/valid;

for k = 1:p.M
    fprintf(' mean clicks (%d)  = %d\n', k,n(k)); 
end
for k = 1:K
    fprintf(' click correlations (%d)  = %d\n', k,Kc(k)); 
end
for k = 1:p.M+1
    fprintf(' probability(%d)  = %d\n', k-1,prob1(k)); 
end

Count = valid;                     % total data records analysed

Nc = n; 
save('expk','Nc','Count');

Kp = Kc;
save('exp_km','Kp','Count');

Cp1 = prob1;
save('expk1','Cp1','Count');

Cp2 = prob2;
save('expk2','Cp2','Count');

Cp4 = prob4;
save('expk4','Cp4','Count');

fprintf(' Total data records analysed = %d\n', Count);
fprintf(' time taken = %d\n', toc()); 
end
