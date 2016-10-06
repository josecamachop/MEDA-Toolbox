function [cumpress,press] = crossval_pca(x,pcs,leave_m,blocks_r,blocks_c,prep,opt)

% Cross-validation for square-prediction-errors computing. The original
% papers are Chemometrics and Intelligent Laboratory Systems 131, 2014, pp.
% 37-50 and Journal of Chemometrics, 26(7), 2012, pp. 361-373.
%
% cumpress = crossval_pca(x,pcs) % minimum call
% [cumpress,press] = crossval_pca(x,pcs,leave_m,blocks_r,blocks_c,prep,opt)
% % complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(x)
%
% leave_m: (str) cross-validation procedure:
%   'rkf': row-wise k fold (default)
%   'ekf': element-wise k fold
%   'cekf': corrected element-wise k fold
%
% blocks_r: [1x1] maximum number of blocks of samples (N by default)
%
% blocks_c: [1x1] maximum number of blocks of variables (M by default)
%
% prep: [1x1] preprocesing
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% opt: (str or num) options for data plotting.
%       0: no plots.
%       1: plot (default)
%
%
% OUTPUTS:
%
% cumpress: [Ax1] Cumulative PRESS
%
% press: [AxM] PRESS per variable.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,8);
% pcs = 0:10;
% cumpress = crossval_pca(X,pcs,'ekf');
%
%
% codified by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Apr/2016
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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
if nargin < 2 || isempty(pcs), pcs = 0:rank(x); end;
A = length(pcs);
if nargin < 3 || isempty(leave_m), leave_m = 'rkf'; end;
if nargin < 4 || isempty(blocks_r), blocks_r = N; end;
if nargin < 5 || isempty(blocks_c), blocks_c = M; end;
if nargin < 6 || isempty(prep), prep = 2; end;
if nargin < 7 || isempty(opt), opt = 1; end;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(leave_m), 'Dimension Error: 3rd argument must be a string. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_c), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
pcs = unique(pcs);

% Validate values of input data
assert (isempty(find(pcs<0)), 'Value Error: 2nd argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(pcs), pcs), 'Value Error: 2nd argumentmust contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (~isempty(strmatch(leave_m,char('rkf','ekf','cekf'))), 'Value Error: 3rd argument must be one of the possible strings. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_r), blocks_r), 'Value Error: 4th argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_c), blocks_c), 'Value Error: 5th argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r>2, 'Value Error: 4th argument must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_c>2, 'Value Error: 5th argument must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r<=N, 'Value Error: 4th argument must be at most N. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_c<=M, 'Value Error: 5th argument must be at most M. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 7th argument must contain a binary value. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Initialization
cumpress = zeros(length(pcs),1);
press = zeros(length(pcs),M);

if ~pcs
    return
end

rows = rand(1,N);
[a,r_ind]=sort(rows);
elem_r=N/blocks_r;

cols = rand(1,M);
[a,c_ind]=sort(cols);
elem_c=M/blocks_c;


% Cross-validation
for i=1:blocks_r,
    
    ind_i = r_ind(round((i-1)*elem_r+1):round(i*elem_r)); % Sample selection
    i2 = ones(N,1);
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
  
    scs = preprocess2Dapp(sample,av,st);
     
    p = pca_pp(ccs,0:max(pcs));
    
    for pc=1:length(pcs),
        
        if pcs(pc) > 0, % PCA Modelling
                               
            p2 = p(:,1:min(pcs(pc),end));
            
            switch lower(leave_m)
                
                case 'rkf',
                    t_est = scs*p2;
                    srec = t_est*p2';
                    pem = sum((scs-srec).^2,1);
                                     
                case 'ekf',
                    t_est = scs*p2;
                    srec = t_est*p2';
                    erec = scs - srec;
                    term3_p = erec;
                    if blocks_c == M,
                        term1_p = (scs-avs_prep).*(ones(ss(1),1)*(sum(p2.*p2,2))');
                    else
                        term1_p = zeros(size(term3_p));
                        for j=1:blocks_c,
                            ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection
                            term1_p(:,ind_j) = (scs(:,ind_j)-avs_prep(:,ind_j))*(p2(ind_j,:)*p2(ind_j,:)');
                        end
                    end
                    
                    term1 = sum(term1_p.^2,1);
                    term2 = sum(2*term1_p.*term3_p,1);
                    term3 = sum(term3_p.^2,1);
                    
                    pem = term1 + term2 + term3;
                    
                case 'cekf'
                    t_cest = ccs*p2;
                    t_sest = scs*p2;
                    
                    rec = t_cest*p2';
                    rec_sam = t_sest*p2';
                    for j=1:blocks_c,
                        ind_j = c_ind(round((j-1)*elem_c+1):round(j*elem_c)); % Variables selection
                        p3 = pca_pp([ccs rec(:,ind_j)],pc);
                        scs2 = [scs rec_sam(:,ind_j)];
                        scs2(:,ind_j) = avs_prep(:,ind_j);
                        t_est = scs2*p3;
                        pred = t_est*p3';
                        srec(:,ind_j) = pred(:,ind_j);
                    end
                    
                    pem = sum((scs-srec).^2,1);
            
                otherwise
                    error('Incorrect leave_m.');
                    
            end
            
        else % Modelling with the average
            pem = sum(scs.^2,1);
        end
        
        press(pc,:) = press(pc,:) + pem;
        
    end            
end

cumpress = sum(press,2);


%% Show results

if opt == '1',    
    fig_h = plot_vec(cumpress,pcs,[],{'#PCs','PRESS'},[],0); 
end

