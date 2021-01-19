function index_fich2 = cfilesys(obslist,centr,label,mult,class,index_fich,thres,path,debug)

% Update of the clustering file system. 
%
% index_fich2 = cfilesys(obslist,centr,label,mult,class,index_fich,thres,path) % minimum call
% index_fich2 = cfilesys(obslist,centr,label,mult,class,index_fich,thres,path,debug) % complete call
%
%
% INPUTS:
%
% obslist: [Nx1] list of groups of observations for the update of the
%   clustering file system. Due to the computation in psc.m, it is assumed
%   that the first observation in a group is the one with highest
%   multiplicity.
%
% centr: [LxM] centroids of the clusters of observations prior to the
%   update.
%
% label: [Lx1] name of the observations (filenames are used by default)
%
% mult: [Lx1] multiplicity of each cluster prior to the update.
%
% class: [Lx1] class associated to each cluster prior to the update.
%
% index_fich: [1xL] cell with the names of the files in the clustering file
%   system prior to the update.
%
% thres: [1x1] maximum number of entries in a file.
%
% path: (str) path to the directory where the output data files are
%   located ('' by default)
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
% index_fich2: [1xN] cell with the names of the files in the clustering file
%   system after the update.
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
assert (nargin >= 8, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(obslist, 1);
[L M] = size(centr);
if nargin < 9 || isempty(debug), debug = 1; end;

% Convert row arrays to column arrays
if size(obslist,1)  == 1, obslist = obslist'; end;
if size(label,1)  == 1, label = label'; end;
if size(mult,1)  == 1, mult = mult'; end;
if size(class,1)  == 1, class = class'; end;

% Convert column arrays to row arrays
if size(index_fich,2)  == 1, index_fich = index_fich'; end;

% Validate dimensions of input data
assert (isequal(size(obslist), [N 1]), 'Dimension Error: 1st argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [L 1]), 'Dimension Error: 3rd argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(mult), [L 1]), 'Dimension Error: 4th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(class), [L 1]), 'Dimension Error: 5th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(index_fich), [1 L]), 'Dimension Error: 6th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(thres), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(debug), [1 1]), 'Dimension Error: 9th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(mult<=0)), 'Value Error: 4th argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(thres<=0)), 'Value Error: 7th argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(thres), thres), 'Value Error: 7th argument must contain an integer. Type ''help %s'' for more info.', routine(1).name);
assert (debug==0 || debug==1 || debug==2, 'Value Error: 9th argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

s = length(obslist);
sc = size(centr,2);

for i=1:s,
    s2 = length(obslist{i});
    indi = obslist{i}(1);
    index_fich2{i} = index_fich{indi};

    if s2 > 1,
        
        indj = find(mult(obslist{i}(2:end))==1);       
        recovered_column = centr(obslist{i}(indj+1),:);
        recovered_label = label(obslist{i}(indj+1));
  
        if ~isempty(recovered_column),
            if mult(indi)>1,
                add_data(index_fich2{i},path,recovered_column,recovered_label,class(indi),'a',thres,[],debug);
            else
                add_data(index_fich2{i},path,[centr(indi,:);recovered_column],{label{indi} recovered_label{:}},class(indi),'w',thres,[],debug);
            end
        end
               
        indj = find(mult(obslist{i}(2:end))>1);
        
        for j=1:length(indj),
            indj2 = obslist{i}(indj(j)+1);
            if mult(indj2)>thres,
                indices = read_indices(index_fich{indj2},path,debug);
                if ispc
                    system(['del ' path index_fich{indj2} '.txt']);
                else
                    system(['rm ' path index_fich{indj2} '.txt']);
                end
                if debug>1, disp(['delete file: ' path index_fich{indj2} '.txt ...']), end;
                add_indices(index_fich2{i},path,indices,debug);
            else
                [recovered_column,recovered_label] = read_data(index_fich{indj2},path,sc,debug);
                if ispc
                    system(['del ' path index_fich{indj2} '.txt']);
                else
                    system(['rm ' path index_fich{indj2} '.txt']);
                end
                if debug>1, disp(['delete file: ' path index_fich{indj2} '.txt ...']), end;
                add_data(index_fich2{i},path,recovered_column,recovered_label,class(indi),'a',thres,[],debug);
            end
        end

    end
end