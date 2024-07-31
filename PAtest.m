function [NoComp, PERCENTILES, Eigs2test]  = PAtest(X,alpha,method,type,corrtype,dispopt)

% Code for performing Horn's Parallel Analysis. Parellel analysis is a well
% establish method to determine the number of Principal components in PCA.
% The performance of this method and its theoretical underpinning have been
% reviwed in  Saccenti, E; Timmerman, M. Considering Horn’s Parallel Analysis 
% from a Random Matrix Theory Point of View Psychometrika 82 (2017)1 
% http://dx.doi.org/10.1007/s11336-016-9515-z. Plese cite this paper if
% using this code
%
% [NoComp, PERCENTILES, Eigs2test]  = PAtest(X,alpha,method,type,corrtype,dispopt)
%       complete call
%
% 
% INPUTS
% 
% X:        n x p data matrix
% alpha:    confidence threshold in ]0 1[
% method:   shuffle:        PA with reshuffling of input data X
%           random:         PA with random data generated under N(0,I)
% type:     covariance:     PA using covariance matrix
%           correlation:    PA using correlation marix
% corrtype  Type of correlation to be used such as spearman or pearson. See
%           Matlab help for correlation type available (type help corr)
% dispopt:  0:              NO graphical display of results
%           >0:              Graphical display of results
% 
%
% OUTPUTS
%
% NoComp:       Number of signifciant component at the alpha elvel scalar)
% PERCENTILES:  (1-alpha) percentiles of the distributions of the 1st,
%               2nd,... eigenvalues of the covariance/correlation matrix
%               (vector)
% Eigs2test:    eigenvalues of the covariance/correlation matrix (vector)
% 
% EXAMPLE OF USE
% 
% SIGMA = eye(10,10);
% SIGMA(1,1) = 5; SIGMA(2,2) = 3;
% MU = zeros(1,10);
% X = mvnrnd(MU,SIGMA,100)
% 
% 
% [NoComp, PERCENTILES, Eigs2test]  = PAtest(X,0.05,'covariance','random','',1)
% [NoComp, PERCENTILES, Eigs2test]  = PAtest(X,0.05,'correlation','random','pearson',1)

% coded by: Edoardo Saccenti (esaccenti@gmail.com)
% last modification: 10/Apr/2024
%
% Copyright (C) 2024  Wageningen Univeristy & Research
% Copyright (C) 2024  Edoardo Saccenti
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

% Main Code

%% Arguments checking

routine=dbstack;
nargin
%assert ((nargin > 6), 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
assert ((nargin == 6) || (nargin == 5) , 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

%% Checking input
if nargin < 6
    dispopt = 0;
end

if alpha > 1 || alpha < 0
    error('Alpha must be in ]0,1[')
end

%% Set some variable
[nobs, nvar] = size(X);
Nrep = 2000;
RAND_EIG = NaN*ones(Nrep, nvar);

%%
disp('Generating random data and percentiles')

switch lower(method)
    
    case 'shuffle' %PA with reshuffling of input data X
        
        for qq = 1 : Nrep
            
            disp(qq)
            
            %Generate "radom" data by shufflig the entries of each column
            %independently. XN columns have same variance of the columns in
            %X.
            XN = NaN*ones(size(X));
            for d = 1:size(X,2)
                XN(:, d) = X(randperm(size(X,1)), d);
            end
            
            
            %Compute covariance/correlation matrix or the "shuffled" data
            switch lower(type)
                case 'covariance'
                    CN = XN'*XN;
                case 'correlation'
                    CN = corr(XN, 'type', corrtype);
                otherwise
                    error('Type not supported')
            end
            
            %Calculate eigenvalues
            RAND_EIG(qq,:) = sort(eig(CN),'descend')';
            
        end
        
    case 'random' %PA with random data generated under N(0,I)
        
        
        for qq = 1 : Nrep
            
            disp(qq)
            
            %Generate radom data with 0 meand Identity covariance N(0,I)
            
            XN = randn(nobs, nvar);
            
            %Compute covariance/correlation matrix
            switch lower(type)
                case 'covariance'
                    CN = XN'*XN;
                case 'correlation'
                    CN = corr(XN, 'type', corrtype);
                otherwise
                    error('Type not supported')
            end
            
            %Calculate eigenvalues
            RAND_EIG(qq,:) = sort(eig(CN),'descend')';
            
        end
        
    otherwise
         error('Method not supported')
        
end

clc
disp('Calculated percentiles')

%% Calculate percentiles for PA
ALPHA = (1-alpha)*100;
PERCENTILES = prctile(RAND_EIG,ALPHA);

%% Compute covariance/correlation of the input data matrix X
switch lower(type)
    case 'covariance'
        C = X'*X;
    case 'correlation'
        C = corr(X, 'type', corrtype);
end

%% Compute sample eigenvalue (eigenvalues of input data matrix X)
Eigs2test = sort(eig(C),'descend')';

%% Print sample eigenvalues and corresponding PA percentiles
for q = 1 : length(Eigs2test)
    
    disp(sprintf('%i. Sample eigenvalue %4.3f - Percentile %4.3f',q,Eigs2test(q),PERCENTILES(q)))
    
end

%% Calculate difference between sample eigenvalues and percentiles and 
%  extract the number of significant components

DeltaEigs = Eigs2test-PERCENTILES;

DeltaEigs(find(DeltaEigs>0)) = 1; 
DeltaEigs(find(DeltaEigs<=0)) = 0;

DeltaEigs(find(DeltaEigs<1,1):end) = 0;

NoComp = sum(DeltaEigs);

disp(sprintf('# of components: %i at alpha %2.2f', NoComp,alpha))

%%
if dispopt > 0 
    %nvar2plot = min(NoComp+10,nvar);
    nvar2plot = nvar;
figure(1)
hold on
grid on; box on
plot(1:nvar2plot,Eigs2test(1:nvar2plot),'.-','MarkerSize',22)
plot(1:nvar2plot,PERCENTILES(1:nvar2plot),'r.-','MarkerSize',22)
xlim([1-0.5 nvar+0.5])
ax = gca;
ax.XTick = [1:1:nvar2plot];
ax.XTickLabel = 1:nvar2plot;
xlabel('Order', 'FontSize',13)
ylabel('Eigenvalue', 'FontSize',13)
legend({sprintf('Original data (%s)',type),sprintf('1-%2.2f percentile %s',alpha, method)},'FontSize',14)

end