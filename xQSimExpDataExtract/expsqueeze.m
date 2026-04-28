function [r] = expsqueeze(p)
%[r] = EXPSQUEEZE(p)
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

%
% Load in Data from available csv files
%

%----------- Start test here to extract data from path -----------

% Get every folder that MATLAB currently has on its search path:
pathList = strsplit(path, pathsep);

% Loop through each entry on the path:
for k = 1:numel(pathList)
    thisFolder = pathList{k};
    if isempty(thisFolder)
        continue
    end

    % Look for *sq.csv files in exactly this folder (no recursion):
    d = dir(fullfile(thisFolder, '*sq.csv'));
    if ~isempty(d)
        fprintf('Squeezing vector found from: %s\n',thisFolder);
        break;
    end
end

%----------- End test here to extract data from path -----------


% d  = dir('*sq.csv');                             %gets files *sq.csv
if isempty(d) 
    error('Found no *sq.csv files in directory\n');
end
if length(d)>1
    error('Found %d *sq.csv files in directory\n', length(d));
end
r = readmatrix(d(1).name);                       %read data from csv file
sz =size(r);
fprintf('Squeezing data is from a %d x %d matrix\n',sz(1),sz(2));
if sz(1) ~= 1
    r = r';
    sz = size(r);
end
if sz(1)~= 1 ||  sz(2) >= p.modes
     error('Largest expected data size is %d x 1 or  1 x %d \n', p.modes,p.modes);
end
if p.verbose ==2
    fprintf('expsqueeze v1_0\n');                %print version name
    for k=1:n
        fprintf('k = %d r = %d\n', k, r(k));     %print data
    end
end
end