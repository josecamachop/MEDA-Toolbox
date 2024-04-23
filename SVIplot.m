function [r2,alpha,q2,res_cross,alpha_cross,betas] = SVIplot(x,varargin)

% Structural and Variance Information plots. The original paper is 
% Chemometrics and Intelligent Laboratory Systems 100, 2010, pp. 48-56. 
%
% r2 = SVIplot(x) % minimum call
% [r2,alpha,q2,res_cross,alpha_cross] = SVIplot(x,'PCs',pcs,'Vars',var,groups,prep,opt) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data setunder analysis
%
% Optional INPUTS (parameter):
%
% 'PCs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(xcs)
%
% 'Vars': (1x1) selected variable for the plot (first variable by default)
%
% 'Groups': [1x1] number of groups in the cross-validation run (7 by default)
%
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% 'Option': (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: no plots
%           1: SVIplot
%       b:
%           0: SVIplot without beta terms
%           1: SVIplot plus beta terms
%   By deafult, opt = '10'. If less than 2 digits are specified, least 
%   significant digit is set to 0, i.e. opt = 1 means a=1 and b=0. If a=0, 
%   then b is ignored.
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
% X = simuleMV(20,10,'LevelCorr',8);
% var = 1;
% [r2,alpha,q2,res_cross,alpha_cross] = SVIplot(X,'PCs',1:3,'Vars',var);
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

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'PCs',0:min(size(x)));   
addParameter(p,'Vars',1);   
addParameter(p,'Groups',7);
addParameter(p,'Preprocessing',2);
addParameter(p,'Option','10');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.PCs;
var = p.Results.Vars;
groups = p.Results.Groups;
opt = p.Results.Option;
prep = p.Results.Preprocessing;


% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
pcs(find(pcs>min(size(x)))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''PCs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(var), [1 1]), 'Dimension Error: parameter ''Vars'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(groups), [1 1]), 'Dimension Error: parameter ''Groups'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==2, 'Dimension Error:parameter ''Option'' must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''PCs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'preprocessing',prep);
p = pca_pp(xcs,'Pcs',1:max(pcs));

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

    [cal_c,m,sd] = preprocess2D(cal,'Preprocessing',prep);
    test_c = preprocess2Dapp(test,m,'SDivideTest',sd);
    
    p = pca_pp(cal_c,'Pcs',1:max(pcs));
    alpha2=0;
    res2 = test_c(:,var);
    for i=1:length(pcs),
        kk = p(:,1:min(pcs(i),size(p,2)));
        q = kk*kk';
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

if opt(1) == '1'
    fig_h=figure;
    hold on
    plot([0 pcs],r2,'.-');
    plot([0 pcs],alpha,'g-x','LineWidth',2);
    plot([0 pcs],q2,'m.-');
    plot(pcs_vect(1:end),alpha_cross(1:end),'ro');
    ch_h=get(fig_h,'Children');
    set(ch_h,'FontSize',14)
    legend('R^2_{A,m}','\alpha^A_{m}','Q^2_{A,m}','\alpha^A_{m}(i)','Location','NorthOutside','Orientation','Horizontal')
    if  opt(2) == '1',
        plot([0 pcs],betas','c')
    end
    
    % Set axis
    axis tight
    ax = axis;
    axis auto
    ax2 = axis;
    axis([ax(1:2) ax2(3:4)])
    
    %legend off
    box on
    hold off

end

