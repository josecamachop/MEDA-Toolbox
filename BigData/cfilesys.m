function index_fich2 = cfilesys(obslist,centr,mult,class,index_fich,thres,path,debug)

% Update of the clustering file system. 
%
% index_fich2 = cfilesys(obslist,centr,mult,class,index_fich,thres,path) % minimum call
% index_fich2 = cfilesys(obslist,centr,mult,class,index_fich,thres,path,debug) % complete call
%
%
% INPUTS:
%
% obslist: (Nx1) list of groups of observations for the update of the
%   clustering file system. Due to the computation in psc.m, it is assumed
%   that the first observation in a group is the one with highest
%   multiplicity.
%
% centr: (LxM) centroids of the clusters of observations prior to the
%   update.
%
% mult: (Lx1) multiplicity of each cluster prior to the update.
%
% class: (Lx1) class associated to each cluster prior to the update.
%
% index_fich: (1xL) cell with the names of the files in the clustering file
%   system prior to the update.
%
% thres: (1x1) maximum number of entries in a file.
%
% path: (str) path to the directory where the output data files are
%   located ('' by default)
%
% debug: (1x1) disply debug messages
%       0: no messages are displayed.
%       1: display only main messages (default) In the present routine, no 
%           messages are displayed. 
%       2: display all messages.
%
%
% OUTPUTS:
%
% index_fich2: (1xN) cell with the names of the files in the clustering file
%   system after the update.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 24/Jan/14.
%
% Copyright (C) 2014  José Camacho Páez
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
    
%% Parameters checking

if nargin < 7, error('Error in the number of arguments.'); end;
if nargin < 8, debug = 1; end;
    
% Computation

s = length(obslist);
sc = size(centr,2);

for i=1:s,
    s2 = length(obslist{i});
    indi = obslist{i}(1);
    index_fich2{i} = index_fich{indi};

    if s2 > 1,
        
        indj = find(mult(obslist{i}(2:end))==1);       
        recovered_column = centr(obslist{i}(indj+1),:);
  
        if ~isempty(recovered_column),
            if mult(indi)>1,
                add_data(index_fich2{i},path,recovered_column,class(indi),'a',thres,[],debug);
            else
                add_data(index_fich2{i},path,[centr(indi,:);recovered_column],class(indi),'w',thres,[],debug);
            end
        end
               
        indj = find(mult(obslist{i}(2:end))>1);
        
        for j=1:length(indj),
            indj2 = obslist{i}(indj(j)+1);
            if mult(indj2)>thres,
                indices = read_indices(index_fich{indj2},path,debug);
                system(['del ' path index_fich{indj2} '.txt']);
                if debug>1, disp(['delete file: ' path index_fich{indj2} ' ...']), end;
                add_indices(index_fich2{i},path,indices,debug);
            else
                recovered_column = read_data(index_fich{indj2},path,sc,debug);
                system(['del ' path index_fich{indj2} '.txt']);
                if debug>1, disp(['delete file: ' path index_fich{indj2} ' ...']), end;
                add_data(index_fich2{i},path,recovered_column,class(indi),'a',thres,[],debug);
            end
        end

    end
end