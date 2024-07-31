function [npcf,rvdim,rvdimperm]=dray(x,varargin)

% Dray's method for PCA component selection. If using this software please cite:
% S. Dray, On the number of principal components: a test of dimensionality
% based on measurements of similarity between matrices, Computational
% Statistics & Data Analysis, 2008 (52), 2228-2237
%
% R. Vitale, J.A. Westerhuis, T. Naes, A.K. Smilde, O.E. de Noord, A.
% Ferrer, Selecting the number of components in Principal Component
% Analysis by permutation testing, Journal of Chemometrics, 2017 (31),
% e2937
%
% npcf = dray(x) % minimum call
% [npcf,rvdim,rvdimperm]=dray(x,'Preprocessing',flagprep,'MaxPermutation',npermmax,'Confidence',conf)
% % complete call
%
%
% INPUTS:
%
% x: [NxM] original data
%
% Optional INPUTS (parameters):
%
% 'Preprocessing': [1x1] 0 for no preprocessing on x; 1 for mean-centering; 2 for
% auto-scaling. Default: 2
%
% 'MaxPermutation': [1x1] maximum number of permutations. Default: 300
%
% 'Confidence': [1x1] confidence level for the test (e.g. 95 for 95%). Default: 95
%
%
% OUTPUTS:
%
% npcf: [1x1] estimated number of PCA factors
%
% rvdim: [1xrank(x)] RVDIM statistic values for all the extractable
% components of x
%
% rvdimperm: [(npermmax)x(npcf+1)] empirical null-distributions of the
% RVDIM statistic for all the significant components plus one.
%
%
% EXAMPLE OF USE: Random data
%
% x = simuleMV(20,10,'LevelCorr',8);
% npcf = dray(x,'MaxPermutation',350,'Confidence',90);
%
%
% codified by: Raffaele Vitale (rvitale86@gmail.com)
% last modification: 8/Apr/2024
%
% Copyright (C) 2024  Universitat Politecnica de Valencia, Valencia
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

%% Parameter checking and initialisation

if nargin<1
    error('Error in the number of arguments. Type ''help dray'' for more info.')
end


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Preprocessing',2);  
addParameter(p,'MaxPermutation',300);
addParameter(p,'Confidence',95);            
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
flagprep = p.Results.Preprocessing;
npermmax = p.Results.MaxPermutation;
conf = p.Results.Confidence;

%% Main code

% Preprocessing

if flagprep==0
    x_p=x;
elseif flagprep==1
    x_p=x-repmat(mean(x),size(x,1),1);
elseif flagprep==2
    x_p=(x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);
end

% Singular Value Decomposition of x

[u,s,v]=svd(x_p);
t=u*s;
eigreal=diag(s.^2)';

rkeff=rank(x_p);
for npc=1:rkeff
    rvdim(npc)=(eigreal(npc))/sqrt(sum(eigreal(npc:rkeff).^2));
end

% Test for the first component

for nperm=1:npermmax
    
    for nvar=1:size(x_p,2)
        xperm(:,nvar)=x_p(randperm(size(x_p,1)),nvar);
    end
    
    [~,sperm,~]=svd(xperm);
    eigperm=diag(sperm.^2);
    rvdimperm(nperm,1)=(eigperm(1))/sqrt(sum(eigperm(1:rkeff).^2));
    
end

% Test for the a-th component (a>1)

while size(rvdimperm,2)<size(rvdim,2) && rvdim(size(rvdimperm,2))>prctile(rvdimperm(:,end),100-(conf/size(rvdimperm,2)))

    e=x_p-t(:,1:size(rvdimperm,2))*v(:,1:size(rvdimperm,2))';
    
    for nperm=1:npermmax
    
        for nvar=1:size(e,2)
            eperm(:,nvar)=e(randperm(size(e,1)),nvar);
        end
        
        [~,sresperm,~]=svd(eperm);
        eigresperm=diag(sresperm.^2);
        rvdimperm_tmp(nperm,1)=eigresperm(1)/sqrt(sum(eigresperm.^2));
        
    end
    
    rvdimperm=[rvdimperm rvdimperm_tmp];
    
end

if size(rvdimperm,2)==size(rvdim,2) && rvdim(size(rvdimperm,2))>prctile(rvdimperm(:,end),100-(conf/size(rvdimperm,2)))
    
    npcf=size(rvdim,2);
    
else

    npcf=size(rvdimperm,2)-1;
    
end

%% Output display

disp(['The number of significant components is ',num2str(npcf)]);