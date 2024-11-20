function model = sparsepls2(X, Y, nlv, nvarX, nvarY, maxiter, tol, mode, mc)
%[max(nlv), unique(nvarX), unique(nvarY), rank(X), rank(Y)]
% Sparse Partial Least Squares. References:
% Lê Cao et al.,Statistical applications in genetics and molecular biology. 2008;7.
% Szymanska et al., Journal of Pharmaceutical and Biomedical Analysis (2016), 27, 170-175
% Szymanska et al., Analytical Chemistry (2015), 87, 869-875.
%
%sparsepls
%
%Summary:
%Function to calculate sparsity PLS models based on the spls function as
%written by Sebastien Dejuan (France), written in R. It uses a procedure
%described in: Le Cao er al., A Sparse PLS for variable selection when
%integrating Omics data (2009).
%
%IO:
%X        = Independent data.
%Y        = Dependent data.
%nlv      = Number of latent variables to model.
%nvarX    = Row vector indicating number of independent variables kept per
%           latent variable modeled.
%nvarY    = Row vector indicating number of dependent variables kept per
%           latent variable modeled.
%maxiter = Maximum number of iterations.
%tol      = Convergence criterion.
%mode     = Use regression (1) or canonical (2) mode for modeling.
%mc       = Option for performing mean-centering (1) or not (0).
%
%model    = Datastructure containing the sparsity PLS model. It has the
%           following fields:
%           P: X loadings.
%           Q: Y loadings.
%           T: X scores.
%           U: Y scores.
%           B: Beta.
%           B0: Intercept values B0.
%           input: Contains all input variables as subfields.
%           mean.X: Mean values for the variables in X.
%           mean.Y: Mean values for the variables in Y.
%           std.X: Standard deviations for the variables in X.
%           std.Y: Standard deviations for the variables in Y.
%           loadings.X: Including amatrix.
%           loadings.Y: Including bmatrix.
%           q: Number of classes.
%           p: Number of variables in the data.
%
%Remarks:
%- Function can perform mean-centering, but only on X.
%- Function handles missing values.
%- Function does not perform scaling.
%
%Author(s)
%- 12/05/2014: Written and adopted by E.Szymanska, Radboud University
%              Nijmegen.
%- 21/05/2014: Edited by G.H. Tinnevelt, Radboud University Nijmegen.
%- 22/02/2016: Edited by T.P.J. Offermans, Radboud University Nijmegen.
%
%Function call example and format:
%[model] = sparsepls2(data, class, 10, repmat(200, 1, 10), repmat(2, 1, 100), 500, 1e-06, 1)
%[model] = sparsepls2(X, Y, nlv, nvarX, nvarY, maxiter, tol, mode, mc)



%VALIDATION OF INPUT

%Check size of X and Y matrices:
[n, p] = size(X);
[r, q] = size(Y);
if length(size(X)) ~= 2 %Works only for matrices with more dimensions.
    error('X must be a numeric matrix');
end
if n ~= r
     error('Unequal number of rows in X and Y');
end

%Check other input parameters:
if nlv <= 0
    error(['Invalid number of latent variables: ' num2str(nlv)]);
end
if nlv > p
    warning('Maximum number of latent variables is set');
    nlv=p;
end
if nlv > max(nvarX)
%    warning('Maximum number of latent variables is set');
    nlv=max(nvarX);
    nvarX(nlv+1:end) = [];
    nvarY(nlv+1:end) = [];
end
if length(nvarX) ~= nlv
    error(['Length of keepX must be equal to: ' num2str(nlv)]);
end
if length(nvarY) ~= nlv
    error(['Length of keepY must be equal to: ' num2str(nlv)]);
end
if any(nvarX > p)
    error(['Each component of keepX must be lower or equal than: ' num2str(p)]);
end
if any(nvarY > q)
    disp (['Each component of keepY must be lower or equal than: ' num2str(q)]);
end   

%Set default input parameters if they are not given:
if nargin < 5 
    maxiter = 500;
    tol = 0.000001;
    mode = 1;
end
if nargin < 3
    nvarX = repmat(100,1,nlv);
    nvarY = repmat(100,1,nlv);
end
    
%Eliminate near-zero variance predictors in X:
[nzvind, X] = nzvElimination(X, n, p);



%INITIALIZATION

%Variable names:
varnamesX = 1:p;
varnamesY = 1:q;
samplenames = 1:n;

%Center:
meanY = mean(Y);
meanX = mean(X);
scaleY = (std(Y));
scaleX = (std(X));
if mc == 1
    %Y = Y - repmat(mean(Y), n, 1);
    X = X - repmat(mean(X), n, 1);
end

%Autoscaling (in original code):
%Y = Y ./ repmat(scaleY, n, 1); 
%X = X ./ repmat(scaleX, n, 1);

%Create intial matrices:
Xtemp = X;
Ytemp = Y;
tmatrix = ones(n, nlv);%Right loading vectors of X'*Y, weights for Y.
umatrix = ones(n, nlv);%Left loading vectors of X'*Y, weights for X.
amatrix = ones(p, nlv);
bmatrix = ones(q, nlv);
cmatrix = ones(p, nlv);%Loading vectors for variables of X.
dmatrix = ones(q, nlv);%Loading vectors for variables of Y in regression mode.
ematrix = ones(q, nlv);%Loading vectors for variables of Y in canonical mode.
nones = ones(1, n);
pones = ones(1, p);
qones = ones(1, q);

%Missing values:
naX = missValues(X);
naY = missValues(Y);
if sum(sum(naX)) > 0
    isnaX = 1;
else
    isnaX = 0;
end
if sum(sum(naY)) > 0
    isnaY = 1;
else
    isnaY = 0;
end



%SPARSE PLS

%Loop for number of components:
for h = 1:nlv
    nx = p - nvarX(:, h);
    ny = q - nvarY(:, h);

    %Step a:
    X1 = Xtemp;
    Y1 = Ytemp;

    %Replace missing values by zeros:
    if isnaX == 1 
        X1(naX == 1) = 0;
    end 
    if isnaY==1 
        Y1(naY == 1) = 0;
    end 
    
    %Step b: svd of crossproduct of X and Y based on simpls:
    M = X1' * Y1;%Crossproduct.
    if find(isnan(M)),
        break;
    end
    [U, ~, V] = svd(M, 0);
    aold = U(:, 1);
    bold = V(:, 1);
    if isnaX == 1 
        t = X1 * aold;
        A = repmat(aold, size(nones));
        A(naX' == 1) = 0;
        anorm = A' * A;
        t = t ./ diag(anorm);
        normt = sqrt(t' * t);
        t = t / normt;
    else
        t = Xtemp * aold / (aold' * aold);
        normt = sqrt(t' * t);
        t = t / normt;
    end
    if isnaY == 1 
        u = Y1 * bold;
        B = repmat(bold, size(nones));
        B(naY' == 1) = 0;
        bnorm = B' * B;
        u = u ./ diag(bnorm);
        normu = sqrt(u' * u);
        u = u / normu;
    else
        u = Ytemp * bold / (bold' * bold);
        normu = sqrt(u' * u);
        u = u / normu;
    end   
    iter = 1;
    
    %Step c: loop for convergence of a and b:
    if isnaX == 1 
        a = X1' * u;
    else 
        a = Xtemp' * u;
    end
    if isnaY == 1 
        b = Y1' * t;
    else 
        b = Ytemp' * t;
    end
    da = a - aold;
    db = b - bold;
 
    while  da' * da > tol %|| db'*db > tol 
        if isnaX == 1 
            a = X1' * u;
        else 
            a = Xtemp' * u;
        end
        if isnaY == 1 
            b = Y1' * t;
        else 
            b = Ytemp' * t;
        end
        if nx ~= 0 
        	[asorted, ~]=sort(abs(a), 1);
            athreshold = asorted(nx, 1);%The smallest coefficient to keep.
        	for z = 1 : size(a, 1)%Number of variables.
            	if abs(a(z, :)) > athreshold
                	tmp = abs(a(z, :)) - athreshold; %Modified GT 21/05/2014.
                	a(z, :) = tmp * sign(a(z, :)); %Keep only variables with coefficients higher than threshold:
                else
                    a(z, :) = 0;
                end
            end
        end
        %a = a / (u' * u);
        norma = sqrt(a' * a);
        a = a / norma;
        if ny ~= 0
        	[bsorted, ~] = sort(abs(b), 1);
            bthreshold = bsorted(ny, 1);%The smallest coefficient to keep.
            for zz = 1:size(b, 1)%Number of variables.
              	if abs(b(zz, :)) > bthreshold
                  	tmp = abs(b(zz, :))- bthreshold;%Keep only variables with coefficients higher than threshold.
                   	b(zz, :) = tmp * sign(b(zz, :)); %Modified GT 21/05/2014.
                else
                    b(zz, :) = 0;
                end
            end
        end
        b = b / (t' * t);
        if isnaX == 1 
            t = X1 * a;
            A = repmat(a, size(nones));
            A(naX' == 1) = 0;
            anorm = A' * A;
            t = t ./ diag(anorm);
            normt = sqrt(t' * t);
            t = t / normt;
        else
            t = Xtemp * a / (a' * a);
            normt = sqrt(t' * t);
            t = t / normt;
        end
        if isnaY == 1 
            u = Y1 * b;
            B = repmat(b, size(nones));
            B(naY' == 1) = 0;
            bnorm = B' * B;
            u = u ./ diag(bnorm);
            normu = sqrt(u' * u);
            u = u / normu;
        else
            u = Ytemp * b / (b' * b);
            normu = sqrt(u' * u);
            u = u / normu;
        end 

        iter = iter + 1;
        da = a - aold;
        db = b - bold;
        aold = a;
        bold = b;
        if iter > maxiter
            %disp('Maximum number of iteration is reached');
            break
        end 
    end
   
    %Step d: deflation of matrices:
    if isnaX == 1 
        X1 = Xtemp;
        X1(naX == 1) = 0;
        c = X1' * t;
        T = repmat(t, size(pones));
        T(naX' == 1) = 0;
        tnorm = T' * T;
        c = c ./ diag(tnorm);
    else
        c = Xtemp' * t / (t' * t);
    end
    Xtemp = Xtemp - (t * c');
    if mode == 1
     if isnaY == 1 %Regression mode.
          Y1 = Ytemp;
        Y1(naY == 1) = 0;
        d = Y1' * t;
        T = repmat(Y1, size(qones));
        T(naY' == 1) = 0;
        tnorm = T' * T;
        d = d ./ diag(tnorm);
     else
         d = Ytemp' * t / (t' * t);
     end
     Ytemp = Ytemp - (t * d');
     e = ematrix(:, h);
    else 
        if isnaY == 1 %Canonical mode.
        Y1 = Ytemp;
        Y1(naY == 1) = 0;
        e = Y1' * u;
        U = repmat(Y1, size(qones));
        U(naY' == 1) = 0;
        unorm = U' * U;
        e = e ./ diag(unorm);
        else
            e = Ytemp' * u / (u' * u);
         d = dmatrix(:, h);
        end
    end
    
     %Collection of results:
     tmatrix(:, h) = t;
     umatrix(:, h) = u;
    
     amatrix(:, h) = a;
     bmatrix(:, h) = b;
     cmatrix(:, h) = c;
     dmatrix(:, h) = d;
     ematrix(:, h) = e;
end
 
%CALCULATE BETA
if mode == 1
    lY = dmatrix;
elseif mode == 2
    lY = ematrix;
end
while isnan(rcond(((cmatrix(:, 1:nlv))' * amatrix(:, 1:nlv))))
    nlv=nlv-1;
end
R = amatrix(:, 1:nlv) / ((cmatrix(:, 1:nlv))' * amatrix(:, 1:nlv));
A = (lY(:, 1:nlv))';

tmatrix(:, nlv+1:end) = [];
umatrix(:, nlv+1:end) = [];

amatrix(:, nlv+1:end) = [];
bmatrix(:, nlv+1:end) = [];
cmatrix(:, nlv+1:end) = [];
dmatrix(:, nlv+1:end) = [];
ematrix(:, nlv+1:end) = [];

B = R * A;
B0 = meanY - meanX * B;






%COLLECTION OF OUTPUT

%Old:
%model.cmatrix = cmatrix; %Loading vectors for variables of X.
%model.dmatrix = dmatrix; %Loading vectors for variables of Y in regression mode.
%model.ematrix = ematrix; %Loading vectors for variables of Y in canonical mode.
%model.input.keepX = nvarX; %Input value for keepX.
%model.input.keepY = nvarY; %Input value for keepY.
%model.input.ncomp = nlv; %Input value for nlv.
%model.input.X = X; %Input value for X.
%model.input.Y = Y; %Input value for Y.
%model.Xmean = meanX; %Mean values for the variables in X.
%model.Ymean = meanY; %Mean values for the variables in Y.
%model.Xscale = scaleX; %Standard deviations for the variables in X.
%model.Yscale = scaleY; %Standard deviations for the variables in Y.
%model.variates.X = tmatrix; %Including tmatrix: right loading vectors of X'*Y, weights for Y.
%model.variates.Y = umatrix; %Including umatrix: left loading vectors of X'*Y, weights for X.
%model.loadings.X = amatrix; %Including amatrix.
%model.loadings.Y = bmatrix; %Including bmatrix.
%model.input.mode = mode; %Input value for mode.



%New, edited by T. Offermans:

model.P = cmatrix; %X loadings.
if mode == 1
    model.Q = dmatrix; %Y loadings in regression mode.
elseif mode == 2
    model.Q = e.matrix; %Y loadings in canonical mode.
end
model.T = tmatrix; %X scores.
model.U = umatrix; %Y Scores.
model.B = B; %Beta.
model.B0 = B0; %Intercept values.
model.R = R; %Weight vector.
model.input.X = X; %Input value for X.
model.input.Y = Y; %Input value for Y.
model.input.nlv = nlv; %Input value for nlv.
model.input.nvarX = nvarX; %Input value for nvarX.
model.input.nvarY = nvarY; %Input value for nvarY.
model.input.maxiter = maxiter; %Input value for maxiter.
model.input.tol = tol; %Input value for tol.
model.input.mode = mode; %Input value for mode.
model.input.mc = mc; %If mean-centering was performed or not.
model.mean.X = meanX; %Mean values for the variables in X.
model.mean.Y = meanY; %Mean values for the variables in Y.
model.std.X = scaleX; %Standard deviations for the variables in X.
model.std.Y = scaleY; %Standard deviations for the variables in Y.

%I do not fully understand what these are:
model.loadings.X = amatrix; %Including amatrix.
model.loadings.Y = bmatrix; %Including bmatrix.

model.q = size(Y, 2);
model.p = size(X, 2);

end




%SUBFUNCTIONS

function[nzvind, Xcorr] = nzvElimination(X, n, p)
% Aim: eliminating near-zero variance predictors in X
% a predictor is classified as near-zero variance if the percentage of unique values in the samples is less than 10% 
thresholdnzv=round(n/10);
nzv=ones(1,p);
for i=1:p
    nzv(1,p)=length(unique(X(:,p)));
end
for i=1:p
    if nzv(1,p)<thresholdnzv
        nzvind(1,p)=1;
    else
        nzvind(1,p)=0;
    end
end
Xcorr=X(:,nzvind == 0);
end

%
function [isnaX] =missValues(data)
% Aim: finding indices of missing values
isnaX=zeros(size(data));
for i=1:size(data,1)
    isnaX(i,:) = isnan(data(i,:));
end
end
  