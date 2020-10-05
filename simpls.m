%SIMPLS Basic SIMPLS.  Performs no error checking.
function [beta,W,P,Q,R] = simpls(X,Y,lvs)

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[N,M] = size(X);
O = size(Y, 2);
if nargin < 3 || isempty(lvs), lvs = 0:rank(X); end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
lvs(find(lvs>M)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (isequal(size(X), [N M]), 'Dimension Error: 1st argument must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Y), [N O]), 'Dimension Error: 2nd argument must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: 3rd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ncomp = max(lvs);
[n,dx] = size(X);
dy = size(Y,2);

% Preallocate outputs
outClass = superiorfloat(X,Y);
Xloadings = zeros(dx,ncomp,outClass);
Yloadings = zeros(dy,ncomp,outClass);
Weights = zeros(dx,ncomp,outClass);


% An orthonormal basis for the span of the X loadings, to make the successive
% deflation X'*Y simple - each new basis vector can be removed from Cov
% separately.
V = zeros(dx,ncomp);

Cov = X'*Y;
for i = 1:ncomp
    % Find unit length ti=X*ri and ui=Y*ci whose covariance, ri'*X'*Y*ci, is
    % jointly maximized, subject to ti'*tj=0 for j=1:(i-1).
    [ri,si,ci] = svd(Cov,'econ'); ri = ri(:,1); ci = ci(:,1); si = si(1);
    ti = X*ri;
    normti = norm(ti); ti = ti ./ normti; % ti'*ti == 1
    Xloadings(:,i) = X'*ti;
    
    qi = si*ci/normti; % = Y'*ti
    Yloadings(:,i) = qi;
    
    Weights(:,i) = ri ./ normti; % rescaled to make ri'*X'*X*ri == ti'*ti == 1

    % Update the orthonormal basis with modified Gram Schmidt (more stable),
    % repeated twice (ditto).
    vi = Xloadings(:,i);
    for repeat = 1:2
        for j = 1:i-1
            vj = V(:,j);
            vi = vi - (vj'*vi)*vj;
        end
    end
    vi = vi ./ norm(vi);
    V(:,i) = vi;

    % Deflate Cov, i.e. project onto the ortho-complement of the X loadings.
    % First remove projections along the current basis vector, then remove any
    % component along previous basis vectors that's crept in as noise from
    % previous deflations.
    Cov = Cov - vi*(vi'*Cov);
    Vi = V(:,1:i);
    Cov = Cov - Vi*(Vi'*Cov);
end

% Postprocessing
R = Weights*inv(Xloadings'*Weights);
W = Weights(:,lvs);
P = Xloadings(:,lvs);
Q = Yloadings(:,lvs);
R = R(:,lvs);
beta=R*Q';