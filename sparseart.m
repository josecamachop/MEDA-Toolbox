
function [AngD, VarA] =  sparseart(X,sparseLoadings)

% Function to calculate the AngD and VarH statistics. The original paper is
% Camacho, Smilde, Saccenti, Westerhuis and Bro. All Sparse PCA Models Are 
% Wrong, But Some Are Useful. Part II: Limitations and Problems of Deflation 
%
% [AngD, VarA] = sparseart(X,sparseLoadings) % complete call
%
%
% INPUTS:
%
% X: [NxM] preprocessed data matrix.
%
% sparseLoadings: [MxA] sparse loading matrix from any sparse PCA alogorithm. 
%       A is the number of components.
%
%
% OUTPUTS:
%           
% AngD: [1xA] The angle (in degrees) of the A sparse loadings to the data 
%       row-space of X. Takes values between 0 and 90 degrees.
%
% VarA:[1xA] Amount of artificially created data (spurious variance not 
%       present in the original data) contaminating components due to the 
%       deflation (null for first component).
%           
%
% EXAMPLE OF USE: Non-overlapping spectra (copy, paste and enjoy)
%
% groups = {1:10,11:20}; % groups for sparse loadings
% vars = 20;
% obs = 5;
% gp = [1 2];
% P = [];
% for i=1:length(gp), % creating sparse loadings
%     pini2 = zeros(vars,1);
%     pini2(groups{gp(i)}) = [.2 .3 .5 .7 .9 .9 .7 .5 .3 .1];
%     P(:,i) = pini2;
%     P(:,i) = P(:,i)/norm(P(:,i));
% end
% 
% for i=1:min(obs, vars),
%     c(i) = 1/(2^i);
% end
% 
% T = repmat(c(1:length(gp)),obs,1).*[[1 1 1 1 0];[1 0 1 0 1]]';
% X = T*P';
% 
% f = figure; subplot(2,1,1), hold on, plot(T,'o-'),
% title('Simulated','FontSize',20), subplot(2,1,2), hold on, bar(P),axis tight 
% 
% [p,t,bel] = gpca(X,groups,'Pcs',1:2);
% 
% [AngD, VarA] = sparseart(X,p)
% 
% f = figure; subplot(2,1,1), hold on, plot(t,'o-'), 
% title('GPCA','FontSize',20), subplot(2,1,2), hold on, bar(p),axis tight 
%
%
% EXAMPLE OF USE: Overlapping spectra (copy, paste and enjoy)
%
% groups = {1:10,6:15,11:20}; % groups for sparse loadings
% vars = 20;
% obs = 5;
% gp = [1 2 3];
% P = [];
% for i=1:length(gp), % creating sparse loadings
%     pini2 = zeros(vars,1);
%     pini2(groups{gp(i)}) = [.2 .3 .5 .7 .9 .9 .7 .5 .3 .1];
%     P(:,i) = pini2;
%     P(:,i) = P(:,i)/norm(P(:,i));
% end
% 
% for i=1:min(obs, vars),
%     c(i) = 1/(2^i);
% end
% 
% T = repmat(c(1:length(gp)),obs,1).*[[1 1 1 1 0];[1 0 1 0 1];[0 1 0 1 0]]';
% X = T*P';
% 
% f = figure; subplot(2,1,1), hold on, plot(T,'o-'),
% title('Simulated','FontSize',20), subplot(2,1,2), hold on, bar(P),axis tight 
% 
% [p,t,bel] = gpca(X,groups,'Pcs',1:3);
% 
% [AngD, VarA] = sparseart(X,p)
% 
% f = figure; subplot(2,1,1), hold on, plot(t,'o-'), 
% title('GPCA','FontSize',20), subplot(2,1,2), hold on, bar(p),axis tight 
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/20024
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

% Get data dimensionality
[obs, vars] = size(X);

% Validate dimensions of input data
assert (isequal(size(sparseLoadings,1), vars), 'Dimension Error: parameter ''sparseLoadings'' must be A-by-M. Type ''help %s'' for more info.', routine(1).name); 

% Get number of components
pcs = 1:size(sparseLoadings,2);

%% Main code

% Deleted non-used columns
ind = find(~sum(abs(sparseLoadings),2));
X(:,ind) = [];
sparseLoadings(ind,:) = [];
vars = size(X,2);

% Perfoma a standard PCA on the data to get standard PCA score T
[~, T] = pca_pp(X);

% Define spase loadings and scores
t = X*sparseLoadings*pinv(sparseLoadings'*sparseLoadings);

% Sign correction
Ploadings = repmat(sign(diag(T'*t))',vars,1).*sparseLoadings;
Tscores = repmat(sign(diag(T'*t))',obs,1).*t;

% Initialize output
VarA = zeros(size(sparseLoadings,2),1);
AngD = zeros(size(sparseLoadings,2),1);

% Calculate AngD and VarA statisrics
for i = 1 : size(sparseLoadings,2)

    a = pinv(X')*Ploadings(:,i);
    AngD(i,1) = round(acosd(dot(Ploadings(:,i),X'*a)/norm(X'*a)));
    
    if i<max(pcs),
        X1 = (X-Tscores(:,1:i)*Ploadings(:,1:i)');
        R = X1*Ploadings(:,i+1)*Ploadings(:,i+1)';
        VarA(i+1,1) = round(100*norm((eye(size(X,2)) - X'*pinv(X'))*(X1'*pinv(X1'))*(R)')/norm(R)); 
    end
    
end

end