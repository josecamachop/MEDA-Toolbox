function omeda_vec = omeda_Lpls(Lmodel,lvs,Ltest,dummy,opt,label)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PCA. The original paper is Journal of Chemometrics, 
% DOI: 10.1002/cem.1405. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda_Lpls(Lmodel,lvs,Ltest,dummy) % minimum call
% omeda_vec = omeda_Lpls(Lmodel,lvs,Ltest,dummy,opt,label) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Ltest.centr: (NxM) centroids of the clusters of observations
%       Ltest.multr: (Nx1) multiplicity of each cluster.
%
% lvs: (1xA) Principal Components considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% Ltest: (struct Lmodel) model with test data:
%       Ltest.XX: (MxM) X-block cross-product matrix.
%       Ltest.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Ltest.centr: (NxM) centroids of the clusters of observations
%       Ltest.multr: (Nx1) multiplicity of each cluster.
%
% dummy: (Nx1) dummy variable containing 1 for the observations in the
%   first group (in test) for the comparison performed in oMEDA, -1 for the 
%   observations in the second group, and 0 for the rest of observations.
%   Also, weights can be introduced.
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
% FALTA CAMBIAR ELLIMITE DE CONTROL DE ACUERDO A omeda_Lpca


%% Parameters checking

if nargin < 4, error('Error in the number of arguments.'); end;
s = size(Lmodel.XX);
if ndims(dummy)==2 & find(size(dummy)==max(size(dummy)))==2, dummy = dummy'; end
if s(1) ~= s(2) || ndims(Lmodel.XX)~=2, error('Error in the dimension of the arguments.'); end;
st = size(Ltest.centr);
if st(2)~=s(2) || size(dummy,1)~=st(1), error('Error in the dimension of the arguments.'); end;
if nargin < 5, opt = 1; end;
if nargin < 6 || isempty(label)
    label=[]; 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end

%% Main code

Lmodel.lv = max(lvs);
[beta,W,P,Q,R] = Lpls(Lmodel);

omeda_vec = omeda(Ltest.centr,Ltest.multr.*dummy,R,P);

%% Show results

if opt == 1,
    plot_vec(omeda_vec,label,'d^2_A');
elseif opt == 2 | opt == 3,
    s = size(Lmodel.centr);
    ov = zeros(100,s(2));
    dummy2 = 2*rand(s(1),1)-1;
    for i=1:100,
        num=randn(s(1),1);
        [kk,ind]=sort(num);
        ov(i,:) = omeda(Lmodel.centr,Lmodel.multr.*dummy2(ind),R,P);
    end
    dev = sqrt(sum(ov.^2)/100);
    
    if opt==2
        plot_vec(omeda_vec,label,'d^2_A',3*[dev;-dev]);
    else       
        idev = find(dev<(1e-2)*max(dev));
        dev(idev)=(1e-2)*max(dev);
        plot_vec(omeda_vec./(3*dev'),label,'d^2_A',[ones(1,s(2));-ones(1,s(2))]);
    end
        
end
    


        