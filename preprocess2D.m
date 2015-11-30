function [xp,average,scale] = preprocess2D(x,prep,weight)

% Preprocess 2-way data.
%
% [xp,average,scale] = preprocess2D(x)          % for mean centering
% [xp,average,scale] = preprocess2D(x,prep,weight)     % complete call
%
% INPUTS:
%
% x: (NxM) Two-way batch data matrix, N(observations) x M(variables)
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)  
%
% weight: (1xM) weight applied after preprocessing scaling. Set to 1 
% by defect.
%
% OUTPUTS:
%
% xp: (NxM) preprocessed data.
%
% average: (1 x M) sample average according to the preprocessing method.
%
% scale: (1 x M) sample scale according to the preprocessing method.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 29/Oct/15.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 2, prep = 1; end;
if (prep<0||prep>2), error('Incorrect value of prep.'); end;
if nargin < 3 || weight==0, 
    weight = ones(1,s(2));
end

% Computation

if s(1)==1 && prep == 2, prep =1; end;

switch prep,
    
    case 1, % mean centering
        
        nanM = isnan(x);
        anM = 1 - nanM;
        x(find(nanM)) = 0;
        average = sum(x,1)./sum(anM,1);    
        scale = ones(1,s(2));
        xp = x - ones(s(1),1)*average;
        xp(find(nanM)) = nan;
        
    case 2, % Trajectory centering and scaling
        
        nanM = isnan(x);
        anM = 1 - nanM;
        x(find(nanM)) = 0;
        average = sum(x,1)./sum(anM,1);
        xc = x - ones(s(1),1)*average; 
        xc(find(nanM)) = 0;
        scale = sqrt(sum(xc.^2,1)./(sum(anM,1)-1));
        ind = find(scale==0);
        scale(ind) = sqrt(ones(1,length(ind))./(2*sum(anM(:,ind),1)-1)); 
        % use 1 by default may reduce detection of anomalous events 
        % what we do is to infer that we need to double the calibration
        % data to find one single element
        xp = xc./(ones(s(1),1)*scale);
        xp(find(nanM)) = nan;
        
    otherwise, % No preprocessing 
        average = zeros(1,s(2));     
        scale = ones(1,s(2)); 
        xp = x;
end

xp = xp.*(ones(s(1),1)*weight);
