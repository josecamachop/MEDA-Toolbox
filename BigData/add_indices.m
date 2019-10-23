function add_indices(name,path,data,debug)

% Add indices to a file in the clustering file system. 
%
% add_indices(name,path,data) % minimum call
% add_indices(name,path,data,debug) % complete call
%
%
% INPUTS:
%
% name: (str) name of the file.
%
% path: (str) path to the directory where the clustering data files are
%   located.
%
% data: [NxM] observations to include in the file.
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
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/May/2017
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(data, 1);
M = size(data, 2);
if nargin < 4 || isempty(debug), debug = 1; end;

% Validate dimensions of input data
assert (isequal(size(debug), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (debug==0 || debug==1 || debug==2, 'Value Error: 4th argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

s=size(data);   
s=max(s);
file=[path name '.txt'];

if debug>1, disp(['add indices in file: ' file ' ...']), end;

fid=fopen(file,'r+');
a=fscanf(fid,'%d',3);
fseek(fid,0,'bof');
str=sprintf('%d %d %d',1,s+a(2),a(3));
str = [str 12*ones(1,10-length(str))];
fprintf(fid,'%s\n',str);
fseek(fid,0,'eof');

for i=1:s,
    fprintf(fid,'%s\n',data{i});
end

fclose(fid);
    
    
    
    