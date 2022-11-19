function [U] = expmatrix(p)
%[U] = EXPMATRIX(p)
%   Gets the experimental transmission data matrix
%   Input transmission data is assumed to come from csv spreadsheets.
%   Real data ends in re.csv, imaginary data ends in im.csv
%   A nonsquare output matrix is padded with extra rows of zeros
%   Input csv data should be in the working directory used by Matlab.
%   Licensed by Peter D Drummond, (2020).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input data
%
if isnumeric(p)  && p <= 0
    U =  'experimental';
    return;
end
%
% Load in Data from available csv files
%

d     = dir('*re.csv');                       % gets files labelled *E.csv
if isempty(d) 
    error('Found 0 real data files in directory\n');
end
if length(d)>1
    error('Found %d real data files in directory\n', length(d));
end
matrix_re = readmatrix(d(1).name);            % Get data from csv file
   d     = dir('*im.csv');                    % gets files labelled *E.csv
if isempty(d) 
    error('Found 0 imaginary data files in directory\n');
end
if length(d)>1
    error('Found %d imaginary data files in directory\n', length(d));
end
matrix_im = readmatrix(d(1).name);           % Get data from csv file

%
% Convert to complex matrix
%

U = matrix_re  + 1i*matrix_im;
sz =size(U);
if sz(1) < p.N ||  sz(2)~= p.M
     error('Data size expected is %d x %d \n',p.N,p.M);
end
fprintf('Transmission data is from a %d x %d matrix\n',sz(1),sz(2));
if sz(1)<sz(2)
    U(sz(1)+1:sz(2),:)=0.0;
end
U=transpose(U);

if p.print ==2
    fprintf('\nexpmatrix v1_0\n\n');                %print version name
    display (U);
end
end