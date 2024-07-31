function [list,var_l] = dat2mat(path,files,output_dir,opt,vars)

% Conversion to dat files output from the parser to mat files which can
% be the input to the MEDA toolbox. We expect that the first column contains
% the names of the variables and the first row the IDs of the observations.
%
% dat2mat(path,files,output_dir,opt,vars) % complete call
%
%
% INPUTS:
%
% path: (str) path where the input files are.
%
% files: (str) regular expresion for the input files.
%
% output_dir: (str) output directory for .mat files.
%
% opt: (1x1) 0 for observations labels identifed as strings and 1 as
%   numbers.
%
% vars: (1xM) variables in the output (all by default)
%
%
% OUTPUTS:
%
% list: {1xL} list with the names of the output files.
%
% var_l: {1xM} list with the names of the variables.
%
% Finally, mat files are stored in the output directory.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Jan/17
%
% Copyright (C) 2024  University of Granada, Granada
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

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, opt = 0; end;
    
	
%% Main code

d = dir(strcat(path,'/',files));

if exist(output_dir) ~= 7 && exist(output_dir) ~= 5,
    mkdir(output_dir);
end
list = {};
for i=1:length(d),
    data = importdata(strcat(path,'/',d(i).name));
    x = data.data;
    
    if opt
        obs_l = num2cell(x(:,1));
        x = x(:,2:end);
    else
        obs_l = data.textdata(2:end,1);
    end
    
    for j=1:length(obs_l)
        obs_l{j} = strtrim(strrep(obs_l{j}, '''', ''));  % remove character '
    end
    
    if i==1,
        s = size(x);
        if nargin < 5, vars = 1:s(2); end
    end
    
    list{end+1} = strcat(output_dir,'/',d(i).name(1:end-4),'.mat');
        
    class = ones(size(x,1),1);
    var_l = data.textdata(1,2:end);
    
    for j=1:length(var_l)
        var_l{j} = strtrim(strrep(var_l{j}, '''', ''));  % remove character ' 
    end
        
    x = x(:,vars);
    var_l = var_l(vars);
    
    disp(strcat('Saving file: ',list{end}))
    save(list{end},'x','class','obs_l','var_l');
end

