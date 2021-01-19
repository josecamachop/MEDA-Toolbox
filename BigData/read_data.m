function [data,label,class,lev,s]=read_data(name,path,nvars,debug)

% Read the data from a file in the clustering file system.  
%
% [data,label,class,lev,s]=read_data(name,path,nvars) % minimum call
% [data,label,class,lev,s]=read_data(name,path,nvars,debug) % complete call
%
%
% INPUTS:
%
% name: (str) name of the file.
%
% path: (str) path to the directory where the clustering data files are
%   located.
%
% nvars: [1x1] number of variables in the data.
%
% debug: [1x1] disply debug messages
%       0: no messages are displayed.
%       1: display only main messages (default) In the present routine, no 
%           messages are displayed.
%       2: display all messages.
%
%
% OUTPUTS:
%
% data: [sxM] observations in the file.
%
% label: [sx1] name of the observations.
%
% class: [1x1] class associated to the observations.
%
% lev: [1x1] hierarchy level of the file, 0 for data and 1 for indices.
%
% s: [1x1] number of observations.
% 
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 12/Jan/2021
%
% Copyright (C) 2021  University of Granada, Granada
% Copyright (C) 2021  Jose Camacho Paez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 4 || isempty(debug), debug = 1; end;

file=[path name '.txt'];
if debug>1, disp(['read data in file: ' file ' ...']), end;

fid=fopen(file,'r');

% Validate dimensions of input data
assert (isequal(size(nvars), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(debug), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isequal(fix(nvars), nvars), 'Value Error: 3rd argument must contain an integer. Type ''help %s'' for more info.', routine(1).name);
assert (debug==0 || debug==1 || debug==2, 'Value Error: 4th argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

a = fscanf(fid,'%d',3);

lev=a(1);
s=a(2);
class=a(3);

if lev==0,
    data = zeros(s,nvars);
    for i=1:s,
        label{i}=fscanf(fid,'%s:',1);
        label{i}(end)=[]; 
        a=fscanf(fid,'%s',1);
        b = textscan(a,'%f',nvars,'Delimiter',',');
        data(i,:) = cell2mat(b);
    end
    fclose(fid);
elseif lev == 1,
    fclose(fid);
    [data,class,lev,s]=read_indices(name,path,debug);
    label = {};
else
    error('Error in the level.');
end
    

