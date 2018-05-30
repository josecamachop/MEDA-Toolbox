function [DImap,best] = DI(T,classes,reg,opt)

% Discriminant Index fro selection of best visualization subspace. The original
% paper is .... 
%
% DImap = DI(T,classes) % minimum call
% [DImap,best] = DI(T,classes,reg,opt) % complete call
%
%
% INPUTS:
%
% T: [NxA] scores of the projection model. 
%
% classes: [Nx1] groups of observations. 0 values will be treated as the
%   rest, and no DImap will be computed for those.
%
% reg: [1x1] regularization parameter, to favour optimums in 1 LV space.
%
% opt: [1x1] options for data plotting:
%       0: no plots
%       1: plot DI matrix per class
%
%
% OUTPUTS:
%
% DImap: [AxAxC] DI matrix per class with code above 0.
%
% best:[Cx2] best subspace for discrimination per class.
%
%
% EXAMPLE OF USE: Seriation and discarding uninformative variables
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 24/May/18
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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
if nargin < 3 || isempty(reg), reg = .1; end; 
if nargin < 4 || isempty(opt), opt = 1; end; 

if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (isequal(size(classes), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(opt<0 | opt>2)), 'Value Error: 3rd argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ucl = unique(classes);
ucl(find(ucl==0))=[];
classesD = zeros(length(classes),length(ucl));
for i=1:length(ucl),
    ind = find(classes==ucl(i));
    classesD(ind,ucl(i)) = 1;
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
                [beta,W,P2,Q,R] = kernel_pls(XX,XY,1);
                T2 = T(:,[i j])*R;
            end
        
            ind = find(classesD(:,k));
            ind2 = find(~classesD(:,k));
            if i==j,
                [kk,m1,dt] = preprocess2D(T(ind,i),2);
                [kk,m2,dt2] = preprocess2D(T(ind2,i),2);
            else
                [kk,m1,dt] = preprocess2D(T2(ind),2);
                [kk,m2,dt2] = preprocess2D(T2(ind2),2);
            end
            
            WS(i,j,k) = dt.^2*(length(ind)-1)+dt2.^2*(length(ind2)-1);
            BS(i,j,k) = (m1 - m2)^2;
       
        end
    end
end

for k=1:length(ucl),    
    DImap(:,:,k) = BS(:,:,k)./WS(:,:,k) .* (ones(size(T,2))+reg*eye(size(T,2))); % I prefer to select single LVs if the difference is below 10%
    
    [topx,topy] = find(DImap(:,:,k)==max(max(squeeze(DImap(:,:,k)))),1);
    best(k,:) = [topx,topy];
end
    
%% Show results

if opt==1,
    for k=1:length(ucl),
        Mv = max(max(squeeze(DImap(:,:,k))));
        plot_map(DImap(:,:,k),[],[0,Mv]);
        ylabel('#LV','FontSize',18)
        xlabel('#LV','FontSize',18)
    end  
end
    
