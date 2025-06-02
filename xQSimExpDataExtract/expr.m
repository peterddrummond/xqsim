function [r] = ExpSqueeze(n)
%[r] = EXPSQUEEZE(N)
%   Gets the experimental squeezing vector r
%   Input squeezing data is assumed to come from csv spreadsheets.
%   Input data should be saved in the working directory used by Matlab
%   Data files should end in sq.csv.
%   Data size should be n x 1 or 1 x n , returns 1 x n squeezing vector.
%   Licensed by Peter D Drummond, (2020).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input data
%

fprintf('ExpSqueeze1_0\n');                 % print current version name

%
% Load in Data from available csv files
%

d     = dir('*sq.csv');                     % gets files labelled *sq.csv
if isempty(d) 
    error('Found 0 squeezing vector files in directory\n');
end
if length(d)>1
    error('Found %d squeezing vector files in directory\n', length(d));
end
r = readmatrix(d(1).name);                  % Get data from csv file
sz =size(r);
fprintf('Data is %d x %d matrix\n',sz(1),sz(2));
if sz(1) ~= 1
    r = r';
    sz = size(r);
end
if sz(1)~= 1 ||  sz(2)~= n
     error('Data size expected is %d x 1 or  1 x %d \n', n,n);
end
for k=1:n
    fprintf('k = %d r = %d\n', k, r(k));
end
end