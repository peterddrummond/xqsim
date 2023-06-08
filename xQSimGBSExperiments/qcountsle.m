function [Count] = qcountsle(p)
%[C] = QCOUNTSLE()
%   Gets the experimental count probability for a GBS experiment
%   Input data must come from little endian binary files
%   Data ends in .bin, last data point is checked for a validity flag
%   Mask is applied to remove time-stamps and specific invalid data points
%   Input data should be saved in the working directory used by Matlab.
%   Produces individual data files for each observable. 
%   Licensed by Peter D. Drummond & Alexander S. Dellios, (2023).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Removes and relabels invalid data, see data description
% 
%

tic()
fprintf('ExpCounts1_0\n');                  % print current version name
maxrecords = 60000000;
%maxrecords = 600000;
fprintf('Processing up to %d records\n', maxrecords); 
%maxrecords = 600;

%
% Mask data for little endian ordering. This section is user specified for
% each data set. 
%

mask = zeros(192,1);
index  = mask;
mask(9:80) = 1;
mask(17) = 0;
mask(36) = 0;
mask(60) = 0;
mask(102:176) = 1;
 
for j=1:192
    index(j) = sum(mask(1:j));
end
index(177:192) = index(177:192) + 1;


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
    elseif sz(2) == 6
        count6 = zeros(p.Part{i}+1);
    end 
end 
click = zeros(p.M,1);
n = zeros(p.M,1);
K = 5;
Kc = zeros(K,1);
K2 = zeros(nchoosek(p.M,2),1);   %All second order correlations
K3 = zeros(nchoosek(p.M,3),1);   %All third order correlations

K2s = zeros(p.M-1,1);           %Subset of second order correlations
K3s = zeros(p.M-2,1);           %Subset of third-order correlations
K4s = zeros(p.M-3,1);           %Subset of fourth-order correlations

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

B2 = cell(1,6);
ln = length(B2);

while bits == p.bits && record < maxrecords 
    for i = 1:ln
        [B2{i},bits] = fread(id,p.bits,'ubit1','ieee-le');
        B2{i} = flipud(B2{i});
    end 
    B1 = cat(1, B2{:});
    if bits == p.bits
        record = record + 1;
        click(1:75)    = B1(176:-1:102);%%75
        click(76:95)   = B1(80:-1:61);%%20
        click(96:118)  = B1(59:-1:37);%%23
        click(119:136) = B1(35:-1:18);%%18
        click(137:144) = B1(16:-1:9);%%8
        n             = n+click;
        
        %All combinations of second order correlations
        A1 = zeros(nchoosek(p.M,2),1);
        s = 1;
        for j=1:p.M-1
            for k=j:p.M-1
                A1(s) = click(j).*click(k+1);
                s = s+1;
            end
        end
        K2 = K2 + A1;
        
        %All combinations of third order correlations
        A2 = zeros(nchoosek(p.M,3),1);
        s = 1;
        for j = 1:p.M-2
            for k = j:p.M-2
                for i = k:p.M-2
                    A2(s) = click(j).*click(k+1,:).*click(i+2,:);
                    s = s+1;
                end 
            end
        end
        K3 = K3 + A2;
        
        
        %Subset of all second order correlations for graphical simplicity
        cln  = length(click);       %Click length
        for i = 1:cln-1
            K2s(i) = K2s(i) + prod(click(i:i+1));
        end 
        
        %Subset of all third order correlations for graphical simplicity
        for i = 1:cln-2
            K3s(i) = K3s(i) + prod(click(i:i+2));
        end 
        
        %Subset of all fourth order correlations for graphical simplicity
        for i = 1:cln-3
            K4s(i) = K4s(i) + prod(click(i:i+3));
        end        
        
        
        
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
        
        c61            = 1+sum(click(1:(1/6*p.M)));
        c62            = 1+sum(click((1/6*p.M+1):(2/6*p.M)));
        c63            = 1+sum(click((2/6*p.M+1):(3/6*p.M)));
        c64            = 1+sum(click((3/6*p.M+1):(4/6*p.M)));
        c65            = 1+sum(click((4/6*p.M+1):(5/6*p.M)));
        c66            = 1+sum(click((5/6*p.M+1):p.M));
        
        count1(c,1) = count1(c,1)+1; 
        count2(c21,c22) = count2(c21,c22)+1;
        count4(c41,c42,c43,c44) = count4(c41,c42,c43,c44)+1;
        count6(c61,c62,c63,c64,c65,c66) = count6(c61,c62,c63,c64,c65,c66) +1;
        valid = valid+1;
    end 
end 

fprintf(' Total records = %d\n',record);
fprintf(' Valid records = %d\n',valid); 

n    = n/valid;
K2s = K2s/valid;
K3s = K3s/valid;
K4s = K4s/valid;

Kc   = Kc/valid;
K2   = K2/valid;
K3   = K3/valid;
prob1 = count1/valid;
prob2 = count2/valid;
prob4 = count4/valid;
prob6 = count6/valid;

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

K2sp = K2s;
save('expk2ms','K2sp','Count');

K3sp = K3s;
save('expk3ms','K3sp','Count');

K4sp = K4s;
save('expk4ms','K4sp','Count');

Kp = Kc;
save('expkm','Kp','Count');

K2p = K2;
save('expk2m','K2p','Count');

K3p = K3;
save('expk3m','K3p','Count');

Cp1 = prob1;
save('expk1','Cp1','Count');

Cp2 = prob2;
save('expk2','Cp2','Count');

Cp4 = prob4;
save('expk4','Cp4','Count');

Cp6 = prob6;
save('expk6','Cp6','Count');

fprintf(' Total data records analysed = %d\n', Count);
fprintf(' time taken = %d\n', toc()); 
end 

