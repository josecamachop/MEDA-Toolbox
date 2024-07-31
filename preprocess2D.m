function [xcs,average,scale] = preprocess2D(x,varargin)

% Preprocess 2-way data.
%
% xcs = preprocess2D(x)          % minimum call
% [xcs,average,scale] = preprocess2D(x,'Preprocessing',prep,'Weight',weights)     % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set
%
% Optional Inputs (Parameters):
% prep: [1x1] preprocesing
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)   
%
% weights: [1xM] weight applied after preprocessing. Set to a vector of 1s 
% by defect.
%
%
% OUTPUTS:
%
% xcs: [NxM] preprocessed data.
%
% average: [1xM] sample average according to the preprocessing method.
%
% scale: [1xM] sample scale according to the preprocessing method.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(10,10,'LevelCorr',8);
% [Xcs,av,sc] = preprocess2D(X);
% fig_h = plot_vec([av' sc'],'XYLabel',{'Avergae','Std Dev'},'Option','11');
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/2024
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Preprocessing',2);   
addParameter(p,'Weight',ones(1,M));
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
prep = p.Results.Preprocessing;
weights = p.Results.Weight;

% Convert column arrays to row arrays
if size(weights,2) == 1, weights = weights'; end;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error in paramaeter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(weights), [1 M]), 'Dimension Error in parameter ''Weights'' be 1-by-M. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (prep>=0 && prep<=3 && isequal(fix(prep), prep), 'Value Error: paramaeter ''Preprocessing'' must contain integers between 0 and 2. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(weights<0)) && isempty(find(weights==Inf)), 'Value Error: parameter ''Weights'' must contain positive values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if N==1 && prep == 2, prep =1; end;

switch prep
    
    case 1 % mean centering
        
        nanM = isnan(x);
        anM = 1 - nanM;
        x(find(nanM)) = 0;
        average = sum(x,1)./sum(anM,1);    
        scale = ones(1,M);
        xcs = x - ones(N,1)*average;
        xcs(find(nanM)) = nan;
        
    case 2 % auto-sclaing
        
        nanM = isnan(x);
        anM = 1 - nanM;
        x(find(nanM)) = 0;
        average = sum(x,1)./sum(anM,1);
        xc = x - ones(N,1)*average; 
        xc(find(nanM)) = 0;
        scale = sqrt(sum(xc.^2,1)./(sum(anM,1)-1));
        ind = find(scale==0);
        scale(ind) = sqrt(ones(1,length(ind))./(2*sum(anM(:,ind),1)-1)); 
        % use 1 by default may reduce detection of anomalous events 
        % what we do is to infer that we need to double the calibration
        % data to find one single element
        xcs = xc./(ones(N,1)*scale);
        xcs(find(nanM)) = nan;
        
    case 3 % pareto-scaling
        
        nanM = isnan(x);
        anM = 1 - nanM;
        x(find(nanM)) = 0;
        average = sum(x,1)./sum(anM,1);
        xc = x - ones(N,1)*average; 
        xc(find(nanM)) = 0;
        scale = sqrt(sqrt(sum(xc.^2,1)./(sum(anM,1)-1)));
        ind = find(scale==0);
        scale(ind) = sqrt(sqrt(ones(1,length(ind))./(2*sum(anM(:,ind),1)-1))); 
        % use 1 by default may reduce detection of anomalous events 
        % what we do is to infer that we need to double the calibration
        % data to find one single element
        xcs = xc./(ones(N,1)*scale);
        xcs(find(nanM)) = nan;
        
    otherwise % No preprocessing 
        average = zeros(1,M);     
        scale = ones(1,M); 
        xcs = x;
end

xcs = xcs.*(ones(N,1)*weights);
