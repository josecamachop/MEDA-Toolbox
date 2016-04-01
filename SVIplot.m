function [r2,alpha,q2,res_cross,alpha_cross] = SVIplot(x,pcs,var,groups,prep,opt)

% Structural and Variance Information plots. The original paper is 
% Chemometrics and Intelligent Laboratory Systems 100, 2010, pp. 48-56. 
%
% r2 = SVIplot(x) % minimum call
% [r2,alpha,q2,res_cross,alpha_cross] = SVIplot(x,pcs,var,groups,prep,opt) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data setunder analysis
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(xcs)
%
% var: (1x1) selected variable for the plot (first variable by default)
%
% groups: [1x1] number of groups in the cross-validation run (7 by default)
%
% prep: [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% opt: [1x1] options for data plotting.
%       0: no plots.
%       1: SVIplot (default)
%       otherwise: SVIplot plus beta terms
%
%
% OUTPUTS:
%
% r2: Goodness of fit.
% 
% alpha: alpha parameter according to the reference.
%
% q2: Goodness of prediction.
%
% res_cross: residuals by CV
%
% alpha_cross: alpha by CV
%
%
% EXAMPLE OF USE: Random data
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% var = 1;
% [r2,alpha,q2,res_cross,alpha_cross] = SVIplot(X,1:3,var);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 29/Mar/2016
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

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 0:rank(x); end;
if nargin < 3 || isempty(var), var = 1; end;
if nargin < 4 || isempty(groups), groups = 7; end; 
if nargin < 5 || isempty(prep), prep = 2; end;
if nargin < 6 || isempty(opt), opt = 1; end; 

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(var), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(groups), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
  
% Validate values of input data
assert (isempty(find(pcs<0)) & isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine.name);
assert (isempty(find(pcs>rank(x))), 'Value Error: 2nd argument must contain values below the rank of the data. Type ''help %s'' for more info.', routine.name);


%% Main code

xcs = preprocess2D(x,prep);
p = pca_pp(xcs,1:max(pcs));

alpha=0;
betas=zeros(M-1,1);
r2 = 0;
for i=1:length(pcs),
    q = p(:,1:pcs(i))*p(:,1:pcs(i))';
    alpha = [alpha q(var,var)];
    betas = [betas q([1:var-1 var+1:end],var)];
    res = xcs*(eye(M)-q);
    r2 = [r2 1-sum(res(:,var).^2)/sum(xcs(:,var).^2)];
end

res_cross=[];
alpha_cross=[];
pcs_vect=[];
rows = rand(1,N);
[a,r_ind]=sort(rows);
elem_r=N/groups;
for j=1:groups,
    ind_i = r_ind(round((j-1)*elem_r+1):round(j*elem_r)); % Sample selection
    i2 = ones(N,1);
    i2(ind_i)=0;
    test = x(ind_i,:);
    cal = x(find(i2),:);
    st = size(test);

    [cal_c,m,sd] = preprocess2D(cal,prep);
    test_c = preprocess2Dapp(test,m,sd);
    
    p = pca_pp(cal_c,1:max(pcs));
    alpha2=0;
    res2 = test_c(:,var);
    for i=1:length(pcs),
        q = p(:,1:pcs(i))*p(:,1:pcs(i))';
        alpha2 = [alpha2 q(var,var)];
        res = test_c*(eye(M)-q);
        res2 = [res2 res(:,var)];
    end
    
    alpha_cross=[alpha_cross;alpha2];
    pcs_vect=[pcs_vect;[0 pcs]];
    res_cross=[res_cross;res2];
    
end

q2 = 1-sum(res_cross.^2)/sum(res_cross(:,1).^2);


%% Show results

if opt,
    fig_h=figure;
    hold on
    plot([0 pcs],r2,'.-');
    plot([0 pcs],alpha,'g-x','LineWidth',2);
    plot([0 pcs],q2,'m.-');
    plot(pcs_vect(1:end),alpha_cross(1:end),'ro');
    ch_h=get(fig_h,'Children');
    set(ch_h,'FontSize',14)
    legend('R^2_{A,m}','\alpha^A_{m}','Q^2_{A,m}','\alpha^A_{m}(i)','Location','NorthOutside','Orientation','Horizontal')
    if opt~=1,
        plot([0 pcs],betas','c')
    end
end

