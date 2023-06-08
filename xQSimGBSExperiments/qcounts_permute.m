function [valid,prob1] = qCOUNTS_permute()
%[valid,prob1] = QCOUNTS()
%   Gets the experimental count probability for a GBS experiment
%   Input data is assumed to come from 128 bit bigendian binary files
%   Data ends in .bin, last data point is checked for a validity flag
%   Mask is applied to remove time-stamps and specific invalid data points
%   Input data should be saved in the working directory used by Matlab.
%   Produces individual data files for each observable.
%   Licensed by Peter D Drummond and Alex Dellios (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Removes and relabels invalid data, see data description
%

tic()
fprintf('qCOUNTS1_0\n');                  % print current version name
maxrecords = 60000000;
%maxrecords = 600000;
fprintf('Processing up to %d records\n', maxrecords); 
%maxrecords = 600;
mask = zeros(128,1);
index = mask;
mask(17:120) = 1;
mask(23) = 0;
mask(25) = 0;
mask(26) = 0;
mask (44) = 0;
mask(128) = 1;

count1 = zeros(101,1);
count2 = zeros(51,51);
count4 = zeros(26,26,26,26);
count5 = zeros(21,21,21,21,21);

click = zeros(100,1);
n = zeros(100,1);
for j=1:120
    index(j) = 101-sum(mask(1:j));
end
K=5;
Kc = zeros(K,1);
M = 100;
K2 = zeros(4950,1);

perm = randperm(100);
%Cp  = cell(7,1);
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
bits = 128;
record = 0;
valid = 0;
while bits  == 128 && record < maxrecords
    [B1,bits] = fread(id,128,'ubit1','ieee-be');
    if bits == 128
      record =record+1;
      if B1(128,1) == 0
        click(1:76)   = B1(120:-1:45);%%76
        click(77:93)  = B1(43:-1:27); %%17
        click(94)     = B1(24);       %%1
        click(95:100) = B1(22:-1:17); %%6
        
        click = click(perm,:);
        
        n             = n+click;
        for i = 1:K
            Kc(i)      = Kc(i)+ prod(click(1:i));
        end
        
        for j=1:M-1
            for k=j:M-1
                A{j}(k-j+1,:) = click(j).*click(k+1);
            end
        end
        K2 = K2 + cat(1,A{:});
        
        c             = 1+sum(click);
        c21            = 1+sum(click(1:50));
        c22            = 1+sum(click(51:100));
        c41            = 1+sum(click(1:25));
        c42            = 1+sum(click(26:50));
        c43            = 1+sum(click(51:75));
        c44            = 1+sum(click(76:100));
        
        c51            = 1+sum(click(1:20));
        c52            = 1+sum(click(21:40));
        c53            = 1+sum(click(41:60));
        c54            = 1+sum(click(61:80));
        c55            = 1+sum(click(81:100));
       
        
        count1(c,1) = count1(c,1)+1; 
        count2(c21,c22) = count2(c21,c22)+1;
        count4(c41,c42,c43,c44) = count4(c41,c42,c43,c44)+1;
        count5(c51,c52,c53,c54,c55) = count5(c51,c52,c53,c54,c55)+1;
        valid = valid+1;
      end
    end
end
fprintf(' Total records = %d\n',record);
fprintf(' Valid records = %d\n',valid); 

n    = n/valid;
Kc   = Kc/valid;
K2   = K2/valid;
prob1 = count1/valid;
prob2 = count2/valid;
prob4 = count4/valid;
prob5 = count5/valid;
for k = 1:100
    fprintf(' mean clicks (%d)  = %d\n', k,n(k)); 
end
for k = 1:K
    fprintf(' click correlations (%d)  = %d\n', k,Kc(k)); 
end
for k = 1:101
    fprintf(' probability(%d)  = %d\n', k-1,prob1(k)); 
end

Count = valid;                     % total data records analysed

Nc = n; 
save('expk','Nc','Count');

Kp = Kc;
save('expkm','Kp','Count');

K2p = K2;
save('expk2m','K2p','Count');

Cp1 = prob1;
save('expk1','Cp1','Count');

Cp2 = prob2;
save('expk2','Cp2','Count');

Cp4 = prob4;
save('expk4','Cp4','Count');

Cp5 = prob5;
save('expk5','Cp5','Count');

Pm = perm;
save('randp','Pm');

fprintf(' Total data records analysed = %d\n', Count);
fprintf(' time taken = %d\n', toc()); 
end