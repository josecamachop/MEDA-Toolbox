function [cumpress,press] = crossval_pls(x,y,lvs,blocks_r,prepx,prepy,opt)

% Row-wise k-fold (rkf) cross-validation for square-prediction-errors computing in PLS.
%
% [cumpress,press] = crossval_pls(x,y,lvs) % minimum call
% [cumpress,press] =
% crossval_pls(x,y,lvs,blocks_r,blocks_c,prepx,prepy,opt) % complete call
%
% INPUTS:
%
% x: (NxM) billinear data set for model fitting
%
% y: (NxO) billinear data set of predicted variables
%
% lvs: (1x1) Latent Variables considered (e.g. pcs = 1:2 selects the
%   first two LVs)
%
% blocks_r: (1x1) maximum number of blocks of samples (Inf by default)
%
% prepx: (1x1) preprocesing of the x-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% prepy: (1x1) preprocesing of the y-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: bar plot (default)
%
% OUTPUTS:
%
% cumpress: (lvs x 1) Cumulative PRESS.
%
% press: (lvs x O) PRESS per variable.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 02/Feb/15.
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


%% Arguments checking

if nargin < 3, error('Error in the number of arguments.'); end;

if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if ndims(y)~=2, error('Incorrect number of dimensions of y.'); end;
sy = size(y);
if find(s<1), error('Incorrect content of y.'); end;

if min(lvs)<0, lvs = 0:max(lvs); end;

if nargin < 4, blocks_r = Inf; end;
if nargin < 5, prepx = 2; end;
if nargin < 6, prepy = 2; end;
if nargin < 7, opt = 1; end;

if blocks_r>s(1), blocks_r = s(1); end
if (blocks_r<2), error('Incorrect value of blocks_r.'); end;


%% Main code

% Initialization
cumpress = zeros(max(lvs)+1,1);
press = zeros(max(lvs)+1,sy(2));

rows = rand(1,s(1));
[a,r_ind]=sort(rows);
elem_r=s(1)/blocks_r;


% Cross-validation
        
for i=1:blocks_r,
    
    ind_i = r_ind(round((i-1)*elem_r+1):round(i*elem_r)); % Sample selection
    i2 = ones(s(1),1);
    i2(ind_i)=0;
    sample = x(ind_i,:);
    calibr = x(find(i2),:); 
    sample_y = y(ind_i,:);
    calibr_y = y(find(i2),:); 

    [ccs,av,st] = preprocess2D(calibr,prepx);
    [ccs_y,av_y,st_y] = preprocess2D(calibr_y,prepy);
        
    scs=sample;
    scs_y=sample_y;
    for j=1:length(ind_i),
        scs(j,:) = (sample(j,:)-av)./st;
        scs_y(j,:) = (sample_y(j,:)-av_y)./st_y;
    end
    
    [beta,W,P,Q,R] = kernel_pls(ccs'*ccs,ccs'*ccs_y,max(lvs));
    
    for lv=lvs,
    
        if lv > 0,
            beta = R(:,max(min(lvs),1):lv)*Q(:,max(min(lvs),1):lv)';
            srec = scs*beta;
            
            pem = scs_y-srec;
            
        else % Modelling with the average
            pem = scs_y;
        end
        
        press(lv+1,:) = press(lv+1,:) + sum(pem.^2,1);
        
    end
end

cumpress = sum(press,2);

%% Show results

if opt == 1,
    fig_h = plot_vec(cumpress(lvs+1),num2str((lvs')),'PRESS',[],1);
    %fig_h = plot_vec2(cumpress(lvs+1),num2str((pcs')),[],'PRESS'); % Problema!!!: plot_vec2 es sólo bar
end

