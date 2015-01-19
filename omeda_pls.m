
function omeda_vec = omeda_pls(cal,y,lvs,test,dummy,prepx,prepy,opt,label,sel)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PLS. The original paper is Journal of Chemometrics, 2011, 25 
% (11): 592-600. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda_pls(cal,y,lvs,test,dummy) % minimum call
% omeda_vec = omeda_pls(cal,y,lvs,test,dummy,prepx,prepy,opt,label) %complete call
%
%
% INPUTS:
%
% cal: (LxM) billinear data set of predictor variables
%
% y: (LxO) billinear data set of predicted variables
%
% lvs: (1xA) Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% test: (NxM) data set with the observations to be compared. These data 
%   should be preprocessed in the same way than calibration data for R and 
%   Q.
%
% dummy: (Nx1) dummy variable containing 1 for the observations in the
%   first group (in test) for the comparison performed in oMEDA, -1 for the 
%   observations in the second group, and 0 for the rest of observations.
%   Also, weights can be introduced.
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
%       1: plot oMEDA vector (default)
%       2: plot oMEDA vector and significance limits
%       3: plot oMEDA vector normalized by significance limits
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
%
% OUTPUTS:
%
% omeda_vec: (Mx1) oMEDA vector.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
%
% Copyright (C) 2014  University of Granada, Granada
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

%% Parameters checking

if nargin < 5, error('Error in the number of arguments.'); end;
s = size(cal);
if ndims(dummy)==2 & find(size(dummy)==max(size(dummy)))==2, dummy = dummy'; end
if s(1) < 1 || s(2) < 1 || ndims(cal)~=2, error('Error in the dimension of the arguments.'); end;
st = size(test);
if st(2)~=s(2) || size(dummy,1)~=st(1), error('Error in the dimension of the arguments.'); end;
if nargin < 6, prepx = 2; end;
if nargin < 7, prepy = 2; end; 
if nargin < 8, opt = 1; end;
if nargin < 9 || isempty(label),
    label=[]; 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end

%% Main code

[cal2,m,sd] = preprocess2D(cal,prepx);
y2 = preprocess2D(y,prepy);

[beta,W,P] = kernel_pls(cal2'*cal2,cal2'*y2,max(lvs));
W2 = W*inv(P'*W);
W2 = W2(:,lvs);
P = P(:,lvs);
        
testp = (test-ones(st(1),1)*m)./(ones(st(1),1)*sd);

omeda_vec = omeda(testp,dummy,W2,P);
    
%% Show results


if opt == 1,
    plot_vec(omeda_vec,label,'d^2_A');
elseif opt == 2 | opt == 3,

    calr = cal2*W2*P';
    l = length(find(dummy));
    vcal = l*sum(calr.^2,1)/s(1);
    lim = 0.9*vcal + 0.1*sum(vcal)/length(vcal); % heuristic
    
    if opt==2
        plot_vec(omeda_vec,label,'d^2_A',[lim;-lim]);
    else       
        omeda_vec = omeda_vec./(lim');
        plot_vec(omeda_vec,label,'d^2_A',[ones(1,s(2));-ones(1,s(2))]);
    end
        
end

        