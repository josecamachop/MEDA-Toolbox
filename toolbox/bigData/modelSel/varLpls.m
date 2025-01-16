function [yvar,tvar] = varLpls(Lmodel,varargin)

% Variability captured in terms of the number of LVs.
%
% varLpls(Lmodel) % minimum call
% [yvar,tvar] = varLpls(Lmodel,'Option',1) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PLS
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.YY: (LxL) Y-block cross-product matrix.
%       Lmodel.lvs: (1x1) number of PCs.
%
% Optional INPUTS (parameters):
%
% 'Plot': (bool) plot results
%       false: no plots.
%       true: plot (default)
%
%
% OUTPUTS:
%
% yvar: [Ax1] Percentage of captured variance of Y.
%
% tvar: [Ax1] Percentage of captured variance of the scores.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Lmodel = iniLmodel(X,Y);
% Lmodel.lvs = 0:10;
% xvar = varLpls(Lmodel);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 16/Jan/2025
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

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Option','1');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
opt = p.Results.Option;

checkLmodel(Lmodel);

% Preprocessing
Lmodel.lvs = unique([0 Lmodel.lvs]);


%% Main code

maxlvs = max(Lmodel.lvs);   
Lmodel = Lpls(Lmodel);
R = Lmodel.altweights;
Q = Lmodel.yloads;

totalVt = sum(eig(Lmodel.XX));
tvar = ones(maxlvs+1,1);
totalVy = sum(eig(Lmodel.YY));
yvar = ones(maxlvs+1,1);
for i=1:maxlvs
    tvar(i+1) = tvar(i+1) - sum(eig(R(:,1:i)'*Lmodel.XX*R(:,1:i)))/totalVt;
    yvar(i+1) = yvar(i+1) - sum(eig(Q(:,1:i)*R(:,1:i)'*Lmodel.XX*R(:,1:i)*Q(:,1:i)'))/totalVy;
end
    
%% Show results

if opt
    plotVec([yvar tvar],'EleLabel',[0 Lmodel.lvs],'XYLabel',{'#LVs','% Residual Variance'},'PlotType','Lines','VecLabel',{'Y','Scores'});
end
