function [data,y,lev,s]=read_indices(name,path,debug)

% Read the data from a file in the clustering file system.  
%
% [data,class,lev,s,indp]=read_data(name,path,nvars) % minimum call
% [data,class,lev,s,indp]=read_data(name,path,nvars,indp,debug) % complete call
%
%
% INPUTS:
%
% name: (str) name of the file.
%
% path: (str) path to the directory where the clustering data files are
%   located.
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
% class: [1x1] class associated to the observations.
%
% lev: [1x1] hierarchy level of the file, 0 for data and 1 for indices.
%
% s: [1x1] number of observations.
% 
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/May/17.
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 3 || isempty(debug), debug = 1; end;

% Validate dimensions of input data
assert (isequal(size(debug), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (debug==0 || debug==1 || debug==2, 'Value Error: 3rd argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

data=[];
file=[path name '.txt'];

if debug>1, disp(['read indices in file: ' file ' ...']), end;

fid=fopen(file,'r');
data={};
cont=1;

a = fscanf(fid,'%d',3);

lev=a(1);
s=a(2);
y=a(3);

for i=1:s,
    data{i} = fscanf(fid,'%s',1);
end   

fclose(fid);

