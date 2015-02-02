function [r2,alpha,q2,res_cross,alpha_cross,fig_h] = SVIplot(data,pcs,var,groups,prep,opt)

% Structural and Variance Information plots. The original paper is 
% Chemometrics and Intelligent Laboratory Systems 100, 2010, pp. 48-56. 
%
% [r2,alpha,q2,res_cross,alpha_cross,fig_h] = SVIplot(data,pcs,var) % minimum call
% [r2,alpha,q2,res_cross,alpha_cross,fig_h] = SVIplot(data,pcs,var,groups,prep,opt) %complete call
%
%
% INPUTS:
%
% data: (NxM) billinear data set under analysis
%
% pcs: (1x1) maximum number of Principal Components considered
%
% var: (1x1) selected variable for the plot.
%
% groups: (1x1) number of groups in the cross-validation run (7 by default)
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: SVIplot (default)
%       2: SVIplot plus beta terms
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
% fig_h: (1x1) figure handle.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
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


if nargin < 3, error('Error in the number of arguments.'); end;
if ndims(data)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(data);
if find(s<1), error('Incorrect content of x.'); end;
if find(var>s(2)), error('Incorrect value of var.'); end;
if nargin < 4, groups = min(s(1),7); end;
if groups < 3, error('Incorrect value of groups.'); end; 
if groups > s(1), groups = s(1); end; 
if nargin < 5, prep = 2; end;
if nargin < 6, opt = 1; end;

data_c = preprocess2D(data,prep);
p = pca_pp(data_c,pcs);

alpha=0;
betas=zeros(s(2)-1,1);
r2 = 0;
for i=1:pcs,
    q = p(:,1:i)*p(:,1:i)';
    alpha = [alpha q(var,var)];
    betas = [betas q([1:var-1 var+1:end],var)];
    res = data_c*(eye(s(2))-q);
    r2 = [r2 1-sum(res(:,var).^2)/sum(data_c(:,var).^2)];
end

if opt,
    fig_h=figure;
    hold on
    plot(0:pcs,r2,'.-');
    plot(0:pcs,alpha,'g-x','LineWidth',2);
end

res_cross=[];
alpha_cross=[];
pcs_vect=[];
rows = rand(1,s(1));
[a,r_ind]=sort(rows);
elem_r=s(1)/groups;
for j=1:groups,
    ind_i = r_ind(round((j-1)*elem_r+1):round(j*elem_r)); % Sample selection
    i2 = ones(s(1),1);
    i2(ind_i)=0;
    test = data(ind_i,:);
    cal = data(find(i2),:);
    st = size(test);

    [cal_c,m,sd] = preprocess2D(cal,prep);
    test_c = (test-ones(st(1),1)*m)./(ones(st(1),1)*sd);
    
    p = pca_pp(cal_c,pcs);
    alpha2=0;
    res2 = test_c(:,var);
    for i=1:pcs,
        q = p(:,1:i)*p(:,1:i)';
        alpha2 = [alpha2 q(var,var)];
        res = test_c*(eye(s(2))-q);
        res2 = [res2 res(:,var)];
    end
    
    alpha_cross=[alpha_cross;alpha2];
    pcs_vect=[pcs_vect;0:pcs];
    res_cross=[res_cross;res2];
    
end

q2 = 1-sum(res_cross.^2)/sum(res_cross(:,1).^2);

if opt,
    plot(0:pcs,q2,'m.-');
    plot(pcs_vect(1:end),alpha_cross(1:end),'ro');
    ch_h=get(fig_h,'Children');
    set(ch_h,'FontSize',14)
    legend('R^2_{A,m}','\alpha^A_{m}','Q^2_{A,m}','\alpha^A_{m}(i)','Location','NorthOutside','Orientation','Horizontal')
    if opt>1,
        plot(0:pcs,betas','c')
    end
end

