function [SDImap,best] = SDI(T,classes,varargin)

% Discriminant Index for selection of best visualization subspace. The original
% paper is Sara Tortorella, Maurizio Servili, Tullia Gallina Toschi, Gabriele 
% Cruciani, José Camacho. Subspace discriminant index to expedite exploration
% of multi-class omics data. Chemometrics and Intelligent Laboratory Systems 206
% (2020) 104160.
%
% SDImap = SDI(T,classes) % minimum call
% [SDImap,best] = SDI(T,classes,'Regular',reg,'Option',opt) % complete call
%
%
% INPUTS:
%
% T: [NxA] scores of the projection model. 
%
% classes: [Nx1] groups of observations. 0 values will be treated as the
%   rest, and no DImap will be computed for those.
%
% Optional INPUTS (parameters):
%
% 'Regular': [1x1] regularization parameter, to favour optimums in 1 LV space
% (0.1 by default)
%
% 'Option': [1x1] options for data plotting:
%       0: no plots
%       1: plot SDI matrix per class (by default)
%
%
% OUTPUTS:
%
% SDImap: [AxAxC] SDI matrix per class with code above 0.
%
% best:[Cx2] best subspace for discrimination per class.
%
%
% EXAMPLE OF USE: Random data with structural relationship  (copy, paste and enjoy)
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% [beta,W,P,Q,R] = simpls(X,Y,'LatVars',lvs);
% T = X*R;
% 
% class = Y;
% class(find(Y==-1))=2;
% SDImap = SDI(T,class);
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

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[N,A] = size(T);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Regular',.1);   
addParameter(p,'Option',1);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
reg = p.Results.Regular;
opt = p.Results.Option;


if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''Classes'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(opt<0 | opt>2)), 'Value Error: parameter ''Option'' must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ucl = unique(classes);
ucl(find(ucl==0))=[];
classesD = zeros(length(classes),length(ucl));
for i=1:length(ucl),
    ind = find(classes==ucl(i));
    classesD(ind,i) = 1;
end

WS = zeros(size(T,2),size(T,2));
BS = zeros(size(T,2),size(T,2));
for i=1:size(T,2),
    for j=1:size(T,2),
        m = zeros(1,length(ucl));
        for k=1:length(ucl),
            
            if i~=j,
                XX=T(:,[i j])'*T(:,[i j]);
                XY=T(:,[i j])'*preprocess2D(classesD(:,k));  
                [beta,W,P2,Q,R] = kernel_pls(XX,XY,'LatVars',1);
                T2 = T(:,[i j])*R;
            end
        
            ind = find(classesD(:,k));
            ind2 = find(~classesD(:,k));
            if i==j,
                [kk,m1,dt] = preprocess2D(T(ind,i),'Preprocessing',2);
                [kk,m2,dt2] = preprocess2D(T(ind2,i),'Preprocessing',2);
            else
                [kk,m1,dt] = preprocess2D(T2(ind),'Preprocessing',2);
                [kk,m2,dt2] = preprocess2D(T2(ind2),'Preprocessing',2);
            end
            
            WS(i,j,k) = dt.^2*(length(ind)-1)+dt2.^2*(length(ind2)-1);
            BS(i,j,k) = (m1 - m2)^2;
       
        end
    end
end

for k=1:length(ucl),    
    SDImap(:,:,k) = BS(:,:,k)./WS(:,:,k) .* (ones(size(T,2))+reg*eye(size(T,2))); % I prefer to select single LVs so I use the regularization parameter
    
    [topx,topy] = find(SDImap(:,:,k)==max(max(squeeze(SDImap(:,:,k)))),1);
    best(k,:) = [topx,topy];
end
    
%% Show results

if opt==1
    for k=1:length(ucl)
        Mv = max(max(squeeze(SDImap(:,:,k))));
       % plot_map(SDImap(:,:,k),[],[0,Mv]);
       plot_map(SDImap(:,:,k),'ColorInt',[0,Mv]);
        ylabel('#LV','FontSize',18)
        xlabel('#LV','FontSize',18)
    end  
end
    
