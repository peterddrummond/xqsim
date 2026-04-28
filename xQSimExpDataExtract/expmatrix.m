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
    d = dir(fullfile(thisFolder, '*re.csv'));
    if ~isempty(d)
        fprintf('Transsmission data real part found from: %s\n',thisFolder);
        break;
    end
end

%----------- End test here to extract data from path -----------


%d     = dir('*re.csv');                       % gets files labelled *E.csv
if isempty(d) 
    error('Found 0 real data files in directory\n');
end
if length(d)>1
    error('Found %d real data files in directory\n', length(d));
end
matrix_re = readmatrix(d(1).name);            % Get data from csv file



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
    d = dir(fullfile(thisFolder, '*im.csv'));
    if ~isempty(d)
        fprintf('Transsmission data imaginary part found from: %s\n',thisFolder);
        break;
    end
end

%----------- End test here to extract data from path -----------



   %d     = dir('*im.csv');                    % gets files labelled *E.csv
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
if sz(1) < length(p.sqz) ||  sz(2)~= p.modes
     error('Data size expected is %d x %d \n',length(p.sqz),p.modes);
end
fprintf('Transmission data is from a %d x %d matrix\n',sz(1),sz(2));
if sz(1)<sz(2)
    U(sz(1)+1:sz(2),:)=0.0;
end
U=transpose(U);

if p.verbose ==2
    fprintf('\nexpmatrix v1_0\n\n');                %print version name
    display (U);
end
end