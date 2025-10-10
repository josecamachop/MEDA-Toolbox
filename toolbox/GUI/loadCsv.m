function vars = loadCsv(filename)

% Load a group of variables from an csv file compatible with the MEDA GUI. 
%
% loadCsv(vars,filename)   % minimum call
%
% See also: saveCsv
%
%
% INPUTS:
%
% filename: (string) name of the Excel file
%
%
% OUTPUTS:
%
% vars: (struct) inputs to store, example
%       vars.X: [NxM] data
%       vars.osb_l: [Nx1] labels of the observations
%       vars.osb_l: [Nx1] labels of the variables
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Get a list of all variables in the workspace
%
% vars = whos;
% 
% % Create an empty struct
% workspace_data = struct();
% 
% % Loop through each variable and add it to the struct
% for i = 1:length(vars)
%     var_name = vars(i).name;
%     % Use dynamic field names to add the variable to the struct
%     workspace_data.(var_name) = eval(var_name);
% end
%
% saveCsv(workspace_data,'MyWorkspace.xlsx');
%
% clear
%
% vars = loadCsv('MyWorkspace.xlsx');
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 23/Aug/2025
% Dependencies: Matlab R2017b, MEDA v1.10
%
% Copyright (C) 2025  University of Granada, Granada
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[~, sheets] = xlsfinfo(filename);

try
    vtype = readmatrix(filename, 'Sheet', 'vtype');
    fprintf('Successfully imported vtype from "%s"\n', filename);
catch ME
    fprintf('Error writing to Excel: %s\n', ME.message);
end

% Pre-allocate a cell array to store the tables
vars = struct();

% Loop through each sheet
for i = 1:numel(sheets)-1
    % Read the data from the current sheet
    sheet_name = sheets{i};

    % Read the data from the sheet
    if vtype(i) == 0
        sheet_table = readmatrix(filename, 'Sheet', sheet_name);
    else
        sheet_table = readcell(filename, 'Sheet', sheet_name);
    end
    
    % Store the table in the structure with the valid name
    vars.(sheet_name) = sheet_table;
    %vars.sheet_name.type

    fprintf('Read and stored sheet: %s\n', sheet_name);
end

end 