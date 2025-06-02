function [Count] = xqcountspnr_permute(p)
%[C] = QCOUNTSPNR()
%   Obtains experimental count probability for a GBS experiment with PNR
%   detectors. 
%   Input data must come from multi-dimensional .mat array of 
%   size n x 1 x M, where M is the number of modes and n is the number of
%   experimental repitions, or shots. 
%   Input data should be saved in the working directory used by Matlab.
%   Produces individual data files for each observable. 
%   Licensed by Peter D Drummond & Alexander S. Dellios, (2023).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic()
fprintf('qCOUNTSPNR1_0\n');
maxrecords = 60000000;
fprintf('Processing up to %d records\n', maxrecords); 


%
% Pre-allocate observables arrays/matrices
%

count1 = zeros(p.max+1,1);                    %Total counts observable
if isfield(p,'Part') == 1                   %Check for partitioned data
    ln = length(p.Part);  
    for i = 1:ln
        szP = size(p.Part{i});
        if szP(2) == 2
            count2 = zeros(p.Part{i}+1);    %2D binning observable
        elseif szP(2) == 4
            count4 = zeros(p.Part{i}+1);    %4D binning observable
        elseif szP(2) == 6
            count6 = zeros(p.Part{i}+1);    %6D binning observable
        end 
    end 
end 
devent = zeros(p.M,1);                      %Detection event, previously 'click'
n = zeros(p.M,1);                           %Mean photon number per channel
n2 = zeros(nchoosek(p.M,2),1);              %All second order correlations

perm = randperm(p.M);         %Generate random permutation

%
% Load in experimental data from available .mat files
%

d     = dir('*pnr.mat');                    %Gets files *pnr.mat
if isempty(d) 
    error('Found 0 experimental data files in directory\n');
end
if length(d)>1
    error('Found %d .mat files\n', length(d));
else
    fprintf('Found %d .mat file\n', length(d)); 
end

es = load(d(1).name);              %Load experimental samples
exp_samples = es.exp_samples;
exp_samples = double(exp_samples);

sz = size(exp_samples);

for i = 1:sz(1)
    s_sample = exp_samples(i,:,:);
    devent = reshape(s_sample, sz(3), 1);

    devent = devent(perm,:);
    
    n = n + devent;
    c = 1 + sum(devent);

    %All combinations of second order correlations
    A1 = zeros(nchoosek(p.M,2),1);
    s = 1;
    for j=1:p.M-1
        for k=j:p.M-1
            A1(s) = devent(j).*devent(k+1);
            s = s+1;
        end
    end
    n2 = n2 + A1;

    if isfield(p,'Part') == 1                   %Check for partitioned data  
        for j = 1:ln
            szP = size(p.Part{j});
            if szP(2) == 2
                c21 = 1+sum(devent(1:(0.5*p.M)));
                c22 = 1+sum(devent((0.5*p.M+1):p.M));
            elseif szP(2) == 4
                c41 = 1+sum(devent(1:(0.25*p.M)));
                c42 = 1+sum(devent((0.25*p.M+1):(0.5*p.M)));
                c43 = 1+sum(devent((0.5*p.M+1):(0.75*p.M)));
                c44 = 1+sum(devent((0.75*p.M+1):p.M));
            elseif szP(2) == 6
                c61 = 1+sum(devent(1:(1/6*p.M)));
                c62 = 1+sum(devent((1/6*p.M+1):(2/6*p.M)));
                c63 = 1+sum(devent((2/6*p.M+1):(3/6*p.M)));
                c64 = 1+sum(devent((3/6*p.M+1):(4/6*p.M)));
                c65 = 1+sum(devent((4/6*p.M+1):(5/6*p.M)));
                c66 = 1+sum(devent((5/6*p.M+1):p.M));
            end 
        end 
    end

    count1(c,1) = count1(c,1)+1; 

    if isfield(p,'Part') == 1                   %Check for partitioned data  
        for j = 1:ln
            szP = size(p.Part{j});
            if szP(2) == 2
                count2(c21,c22) = count2(c21,c22)+1;    
            elseif szP(2) == 4
                count4(c41,c42,c43,c44) = count4(c41,c42,c43,c44)+1;    
            elseif szP(2) == 6
                count6(c61,c62,c63,c64,c65,c66) = count6(c61,c62,c63,c64,c65,c66) +1;    
            end 
        end 
    end 
end 

n    = n/sz(1);
n2   = n2/sz(1);
prob1 = count1/sz(1);
if isfield(p,'Part') == 1               %Check for partitioned data  
    for j = 1:ln
        szP = size(p.Part{j});
        if szP(2) == 2
            prob2 = count2/sz(1);    
        elseif szP(2) == 4
            prob4 = count4/sz(1);    
        elseif szP(2) == 6
            prob6 = count6/sz(1);    
        end 
    end 
end


Count = sz(1);                         %Total data records analysed

Nc = n; 
save('expn','Nc','Count');

n2p = n2;
save('expn2m','n2p','Count');

Cp1 = prob1;
save('expn1','Cp1','Count');

Pm = perm;
save('randp','Pm');

if isfield(p,'Part') == 1               %Check for partitioned data  
    for j = 1:ln
        szP = size(p.Part{j});
        if szP(2) == 2
            Cp2 = prob2;
            save('expn2','Cp2','Count');    
        elseif szP(2) == 4
            Cp4 = prob4;
            save('expn4','Cp4','Count');    
        elseif szP(2) == 6
            Cp6 = prob6;
            save('expn6','Cp6','Count');    
        end 
    end 
end

fprintf(' Total data records analysed = %d\n', Count);
fprintf(' Time taken = %d\n', toc()); 



