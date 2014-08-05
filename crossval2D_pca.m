function [cumpress,press,pem] = crossval2D_pca(x,pc,leave_m,blocks_r,blocks_c,prep)

% Cross-validation for square-prediction-errors computing. The original
% papers are Chemometrics and Intelligent Laboratory Systems 131, 2014, pp.
% 37-50 and Journal of Chemometrics, 26(7), 2012, pp. 361-373.
%
% [cumpress,press,pem] = crossval2D_pca(x,pc) % minimum call
% [cumpress,press,pem] = crossval2D_pca(x,pc,leave_m,blocks_r,blocks_c,scaling)
% % complete call
%
% INPUTS:
%
% x: (NxM) billinear data set for model fitting
%
% pc: (1x1) number of principal components
%
% leave_m: (str) cross-validation procedure:
%   'rkf': row-wise k fold (default)
%   'ekf': element-wise k fold
%   'eekf': same as 'ekf' but in a fast way
%   'eekf2': same as 'ekf' but in a fast way
%   'cekf': corrected element-wise k fold
%
% blocks_r: (1x1) maximum number of blocks of samples (Inf by default)
%
% blocks_c: (1x1) maximum number of blocks of variables (Inf by default)
%
% prep: (1x1) preprocesing
%       0: no preprocessing (default)
%       1: mean-centering 
%       2: auto-scaling 
%
%
% OUTPUTS:
%
% cumpress: (1x1) Cumulative PRESS.
%
% press: (1xM) PRESS per variable.
%
% pem: (NxM) Matrix containing prediction errors.
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
%
% Copyright (C) 2014  José Camacho Páez
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

% Arguments checking

if nargin < 2, error('Error in the number of arguments.'); end;

if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;

if nargin < 3, leave_m = 'rkf'; end;
if nargin < 4, blocks_r = Inf; end;
if nargin < 5, blocks_c = Inf; end;
if nargin < 6, prep = 2; end;

if pc<0, error('Incorrect value of pc.'); end;
if blocks_r>s(1), blocks_r = s(1); end
if (blocks_r<2), error('Incorrect value of blocks_r.'); end;
if blocks_c>s(2), blocks_c = s(2); end
if (blocks_c<2), error('Incorrect value of blocks_c.'); end;


% Initialization

pem = zeros(s);
cumpress = 0;
press = zeros(1,s(2));

if isequal(lower(leave_m), 'cekf'),
    xcs = preprocess2D(x,prep);
    [p,t] = pca_pp(xcs,pc);
end

rows = rand(1,s(1));
[a,r_ind]=sort(rows);
elem_r=s(1)/blocks_r;

cols = rand(1,s(2));
[a,c_ind]=sort(cols);
elem_c=s(2)/blocks_c;


% Cross-validation
for i=1:blocks_r,
    
    ind_i = r_ind(round((i-1)*elem_r+1):round(i*elem_r)); % Sample selection
    i2 = ones(s(1),1);
    i2(ind_i)=0;
    sample = x(ind_i,:);
    calibr = x(find(i2),:); 
    sc = size(calibr);
    ss = size(sample);

    [ccs,av,st] = preprocess2D(calibr,prep);
    
    if ~prep,
        avs_prep=ones(ss(1),1)*mean(ccs);
    else
        avs_prep=zeros(ss);
    end
        
    scs=sample;
    for j=1:length(ind_i),
        scs(j,:) = (sample(j,:)-av)./st;
    end
    srec = zeros(size(scs));
    
    if pc > 0, % PCA Modelling
            
        switch lower(leave_m)
    
            case 'rkf',
                p = pca_pp(ccs,pc);
                t_est = scs*p;
                srec = t_est*p';                    
                
            case 'eekf2',
                p = pca_pp(ccs,pc);
                t_est = scs*p;
                srec = t_est*p';
                pem(ind_i,:) = (scs-avs_prep).*(sum(p.*p,2))' + scs - srec;
                
                
            case 'eekf',
                p = pca_pp(ccs,pc);
                t_est = scs*p;
                srec = t_est*p';
                erec = scs - srec;
                for j=1:blocks_c,
                    ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection                                  
                    pem(ind_i,ind_j) = (scs(:,ind_j)-avs_prep(:,ind_j))*(p(ind_j,:)*p(ind_j,:)') + erec(:,ind_j);
                end
   
            case 'ekf',
                p = pca_pp(ccs,pc);
                for j=1:blocks_c,
                    ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection
                    scs2 = scs;
                    scs2(:,ind_j) = avs_prep(:,ind_j);
                    t_est = scs2*p;
                    pred = t_est*p';
                    srec(:,ind_j) = pred(:,ind_j);                                  
                end

            case 'cekf',             
                [rec,av,st] = preprocess2D(t(find(i2),:)*p',prep);             
                rec_sam=t(ind_i,:)*p';
                for j=1:length(ind_i),
                    rec_sam(j,:) = (rec_sam(j,:)-av)./st;
                end
                for j=1:blocks_c,
                    ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection
                    p2 = pca_pp([ccs rec(:,ind_j)],pc);
                    scs2 = [scs rec_sam(:,ind_j)];
                    scs2(:,ind_j) = avs_prep(:,ind_j);
                    t_est = scs2*p2;
                    pred = t_est*p2';
                    srec(:,ind_j) = pred(:,ind_j);                                  
                end
                
           otherwise
               error('Incorrect leave_m.');

        end
        
        if isequal(lower(leave_m),'ekf') | isequal(lower(leave_m),'rkf') | isequal(lower(leave_m),'cekf'),
            pem(ind_i,:) = scs-srec;
        end
        
    else % Modelling with the average
        pem(ind_i,:) = scs;
    end
    
end

press = sum(pem.^2,1);

cumpress = sum(press);



