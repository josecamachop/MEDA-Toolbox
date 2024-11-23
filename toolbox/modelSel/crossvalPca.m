function [cumpress,press] = crossvalPca(x,pcs,varargin)

% Cross-validation for square-prediction-errors computing. The original
% papers are Chemometrics and Intelligent Laboratory Systems 131, 2014, pp.
% 37-50 and Journal of Chemometrics, 26(7), 2012, pp. 361-373.
%
% cumpress = crossvalPca(x,pcs) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(x)
%
%
% Optional INPUTS (parameters):
%
% 'ValProcedure': (str) cross-validation procedure:
%   'rkf': row-wise k fold (default)
%   'ekf': element-wise k fold
%   'cekf': corrected element-wise k fold
%
% 'MaxSampleBlock': [1x1] maximum number of blocks of samples (N by default)
%
% 'MaxVarBlock': [1x1] maximum number of blocks of variables (M by default)
%
% 'Preprocesing': [1x1] preprocesing
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% 'Option': (bool) plot results.
%       false: no plots.
%       true: plot (default)
%
%
% OUTPUTS:
%
% cumpress: [Ax1] Cumulative PRESS
%
% press: [AxM] PRESS per variable.
%
%
% EXAMPLE OF USE: Random data using mean centering
% 
% X = simuleMV(20,10,'LevelCorr',8);
% pcs = 0:10;
% cumpress = crossvalPca(X,pcs,'ValProcedure','rkf','Preprocessing', 1);
%
%
% codified by: Jose Camacho (josecamacho@ugr.es)
% last modification: 20/Nov/2024
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
% if nargin < 2 || isempty(pcs), pcs = 0:rank(x); end;
A = length(pcs);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'ValProcedure','rkf'); 
addParameter(p,'MaxSampleBlock',N);
addParameter(p,'MaxVarBlock',M);
addParameter(p,'Preprocessing',2);   
addParameter(p,'Plot',true);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
leavem = p.Results.ValProcedure;
blocksr = p.Results.MaxSampleBlock;
blocksc = p.Results.MaxVarBlock;
prep = p.Results.Preprocessing;
opt = p.Results.Plot;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(leavem), 'Dimension Error: parameter ''ValProcedure'' must be a string. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxSampleBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksc), [1 1]), 'Dimension Error: parameter ''MaxVarBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: parameter ''Plot'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
pcs = unique(pcs);

% Validate values of input data
assert (isempty(find(pcs<0)), 'Value Error: parameter ''pcs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(pcs), pcs), 'Value Error: parameter ''pcs'' contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (~isempty(strmatch(leavem,char('rkf','ekf','cekf'))), 'Value Error: parameter ''ValProcedure'' must be one of the possible strings. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxSampleBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksc), blocksc), 'Value Error: parameter ''MaxVarBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>2, 'Value Error: parameter ''MaxSampleBlock'' must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocksc>2, 'Value Error: parameter ''MaxVarBlock'' must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxSampleBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);
assert (blocksc<=M, 'Value Error: parameter ''MaxVarBlock'' must be at most M. Type ''help %s'' for more info.', routine(1).name);
assert (islogical(opt), 'Value Error: parameter ''Plot'' must contain a boolean. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Initialization
cumpress = zeros(length(pcs),1);
press = zeros(length(pcs),M);

if ~pcs
    return
end

rows = rand(1,N);
[a,rind]=sort(rows);
elemr=N/blocksr;

cols = rand(1,M);
[a,cind]=sort(cols);
elemc=M/blocksc;


% Cross-validation
for i=1:blocksr
    
    indi = rind(round((i-1)*elemr+1):round(i*elemr)); % Sample selection
    i2 = ones(N,1);
    i2(indi)=0;
    sample = x(indi,:);
    calibr = x(find(i2),:); 
    sc = size(calibr);
    ss = size(sample);

    [ccs,av,st] = preprocess2D(calibr,'Preprocessing',prep);
    
    if ~prep
        avsprep=ones(ss(1),1)*mean(ccs);
    else
        avsprep=zeros(ss);
    end
  
    scs = preprocess2Dapp(sample,av,'Scale',st);
     
    model = pcaEig(ccs,'PCs',0:max(pcs));
    p = model.loads;
    
    for pc=1:length(pcs)
        
        if pcs(pc) > 0 % PCA Modelling
                               
            p2 = p(:,1:min(pcs(pc),end));
            
            switch lower(leavem)
                
                case 'rkf'
                    test = scs*p2;
                    srec = test*p2';
                    pem = sum((scs-srec).^2,1);
                                     
                case 'ekf'
                    test = scs*p2;
                    srec = test*p2';
                    erec = scs - srec;
                    term3p = erec;
                    if blocksc == M
                        term1p = (scs-avsprep).*(ones(ss(1),1)*(sum(p2.*p2,2))');
                    else
                        term1p = zeros(size(term3p));
                        for j=1:blocksc
                            indj = cind(round((j-1)*elemc+1):round(j*elemc)); % Variables selection
                            term1p(:,indj) = (scs(:,indj)-avsprep(:,indj))*(p2(indj,:)*p2(indj,:)');
                        end
                    end
                    
                    term1 = sum(term1p.^2,1);
                    term2 = sum(2*term1p.*term3p,1);
                    term3 = sum(term3p.^2,1);
                    
                    pem = term1 + term2 + term3;
                    
                case 'cekf'
                    tcest = ccs*p2;
                    tsest = scs*p2;
                    
                    rec = tcest*p2';
                    recsam = tsest*p2';
                    for j=1:blocksc
                        indj = cind(round((j-1)*elemc+1):round(j*elemc)); % Variables selection
                        model3 = pcaEig([ccs rec(:,indj)],pc);
                        p3 = model3.loads;
                        scs2 = [scs recsam(:,indj)];
                        scs2(:,indj) = avsprep(:,indj);
                        test = scs2*p3;
                        pred = test*p3';
                        srec(:,indj) = pred(:,indj);
                    end
                    
                    pem = sum((scs-srec).^2,1);
            
                otherwise
                    error('Incorrect leavem.');
                    
            end
            
        else % Modelling with the average
            pem = sum(scs.^2,1);
        end
        
        press(pc,:) = press(pc,:) + pem;
        
    end            
end

cumpress = sum(press,2);


%% Show results

if opt    
    figh = plotVec(cumpress,'EleLabel',pcs,'XYLabel',{'#PCs','PRESS'},'Option','01'); 
end

