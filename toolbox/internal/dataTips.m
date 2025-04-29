function output_txt = dataTips(~, event_obj, bdata, varargin)
% Data tips for scatter and vec plots
%
% See also: plotScatter, plotVec
%
% INPUTS:
%
% bdata: Data used in the plots.
%   (Nx2) for scatter plots
%   (Nx1) for vec plots
%
% Optional INPUTS (parameters):
%
% 'EleLabel': [Nx1] name of the elements
%
% 'ObsClass': [Nx1] groups for different visualization
%
% 'ClassType': str. Class type of the data (Categorial or numerical). Default = "Categorical".
%
% OUTPUTS:
%
% output_txt: dataTip for point/bar.
%
% coded by: Jesús García Sánchez (gsus@ugr.es)
%           
% last modification: 03/Feb/2025
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

%

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

N = size(bdata, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',[]);   
addParameter(p,'Classes',[]);
addParameter(p,'ClassType',"Categorical");
addParameter(p,'Multiplicity',ones(N,1));
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
classes = p.Results.Classes;
classType = p.Results.ClassType;
mult = p.Results.Multiplicity;

if iscell(classes)
    classes = string(cell2mat(classes)); end

if iscell(elabel)
    elabel = string(cell2mat(elabel)); end

%% Main code

    pos = get(event_obj, 'Position');
    if size(bdata,2) == 1 % Vector
        idx = pos(1);
    end
    if size(bdata,2) == 2 % Scatter
        idx = find(bdata(:,1) == pos(1) & bdata(:,2) == pos(2));
    end
    output_txt = {
        ['x: ', num2str(pos(1))], ...
        ['y: ', num2str(pos(2))],    };

    
    if ~isempty(elabel) || ~isempty(classes) % Newline
        output_txt{end+1} = ''; end

    if ~isempty(elabel)
        output_txt{end+1} = ['Observation: ', num2str(elabel(idx))];
    end

    if ~isempty(classes)
        if classType == "Categorical"
        output_txt{end+1} = ['Class: ', num2str(classes(idx))]; end
        if classType == "Numerical"
        output_txt{end+1} = ['Value: ', num2str(classes(idx))]; end
    end

    if mult ~= ones(N,1)
        output_txt{end+1} = ['Multiplicity: ', num2str(mult(idx))];
    end
end