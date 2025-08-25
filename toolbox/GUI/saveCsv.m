function saveCsv(vars,filename)

% Save a group of variables in an csv file compatible with the MEDA GUI. 
%
% saveCsv(vars,filename)   % minimum call
%
% See also: loadCsv
%
%
% INPUTS:
%
% vars: (struct) inputs to store, example
%       vars.X: [NxM] data
%       vars.osb_l: [Nx1] labels of the observations
%       vars.osb_l: [Nx1] labels of the variables
%
% filename: (string) name of the Excel file
%
%
% OUTPUTS:
%
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

%% Main code

% Check if the file exists
if exist(filename, 'file') == 2
    % If the file exists, delete it
    delete(filename);
end

% Get all field names of the struct
fieldNames = fieldnames(vars);

% Loop through each field name
for i = 1:length(fieldNames)
    % Get the current field name
    fieldName = fieldNames{i};

    % Get the value of that field using dynamic field notation
    fieldValue = vars.(fieldName);

    try
        if iscell(fieldValue)
            vtype(i) = 1;
            writecell(fieldValue, filename, 'Sheet', fieldName);
        else
            vtype(i) = 0;
            writematrix(fieldValue, filename, 'Sheet', fieldName);
        end
        fprintf('Successfully exported %s to "%s"\n', fieldName, filename);
    catch ME
        fprintf('Error writing to Excel: %s\n', ME.message);
    end
end

try
    writematrix(vtype, filename, 'Sheet', 'vtype');
    fprintf('Successfully exported vtype to "%s"\n', filename);
catch ME
    fprintf('Error writing to Excel: %s\n', ME.message);
end

end