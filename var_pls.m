
function [y_var,t_var] = var_pls(x,y,varargin)

% Variability captured in terms of the number of LVs.
%
% var_pls(x,y) % minimum call
% var_pls(x,y,'LatVars',maxlvs,'PreprocessingX',prepx,'PreprocessingY',prepy,'Option',opt) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% Optional INPUTS (parameter):
%
% 'LatVars': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x). The value for 0 PCs is
%   added at the begining if not specified.
%
% 'PreprocessingX': [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% 'PreprocessingY': [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% 'Option': (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: no plots
%           1: plot Residual variance
%       b:
%           0: Residual Variance in Y and Scores 
%           1: Residual Variance in Y 
%   By deafult, opt = '10'. If less than 2 digits are specified, least 
%   significant digit is set to 0, i.e. opt = 1 means a=1 and b=0. If a=0, 
%   then b is ignored.
%
%
% OUTPUTS:
%
% y_var: [Ax1] Percentage of captured variance of Y.
%
% t_var: [Ax1] Percentage of captured variance of the scores.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 0:10;
% x_var = var_pls(X,Y,'LatVars',lvs);
%
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'LatVars',0:rank(x));   
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);   
addParameter(p,'Option',10);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LatVars;
opt = p.Results.Option;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''LatVars'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LatVars'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==2, 'Dimension Error: parameter ''Option'' must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique([0 lvs]);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LatVars'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LatVars'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prepx); 
ycs = preprocess2D(y,'Preprocessing',prepy); 

%[beta,W,P,Q,R] = kernel_pls(xcs'*xcs,xcs'*ycs,1:max(lvs));
[beta,W,P,Q,R] = simpls(xcs,ycs,'LatVars',1:max(lvs));
lvs(find(lvs>size(W,2))) = [];

totalVt = sum(sum(xcs.^2));
t_var = ones(length(lvs),1);
totalVy = sum(sum(ycs.^2));
y_var = ones(length(lvs),1);
for i=1:length(lvs),
    t_var(i) = t_var(i) - sum(eig(R(:,1:lvs(i))'*xcs'*xcs*R(:,1:lvs(i))))/totalVt;
    y_var(i) = y_var(i) - sum(eig(Q(:,1:lvs(i))*R(:,1:lvs(i))'*xcs'*xcs*R(:,1:lvs(i))*Q(:,1:lvs(i))'))/totalVy;
end
    
%% Show results
           
if opt(1) == '1'
    if opt(2) == '1'
        plot_vec(y_var,'EleLabel',lvs,'XYLabel',{'#PCs','% Residual Variance in Y'},'Option','01');
    else
        plot_vec([y_var t_var],'EleLabel',lvs,'XYLabel',{'#PCs','% Residual Variance'},'Option','01','VecLabel',{'Y','Scores'});
        legend('show');
    end
end
