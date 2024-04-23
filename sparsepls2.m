function [model]=sparsepls2(X, Y, nlv, nvarX, nvarY, max_iter, tol, mode, mc)
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
%max_iter = Maximum number of iterations.
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
%           loadings.X: Including a_matrix.
%           loadings.Y: Including b_matrix.
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
%[model] = sparsepls2(X, Y, nlv, nvarX, nvarY, max_iter, tol, mode, mc)



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
    max_iter = 500;
    tol = 0.000001;
    mode = 1;
end
if nargin < 3
    nvarX = repmat(100,1,nlv);
    nvarY = repmat(100,1,nlv);
end
    
%Eliminate near-zero variance predictors in X:
[nzv_ind, X] = nzv_elimination(X, n, p);



%INITIALIZATION

%Variable names:
var_names_X = 1:p;
var_names_Y = 1:q;
sample_names = 1:n;

%Center:
mean_Y = mean(Y);
mean_X = mean(X);
scale_Y = (std(Y));
scale_X = (std(X));
if mc == 1
    %Y = Y - repmat(mean(Y), n, 1);
    X = X - repmat(mean(X), n, 1);
end

%Autoscaling (in original code):
%Y = Y ./ repmat(scale_Y, n, 1); 
%X = X ./ repmat(scale_X, n, 1);

%Create intial matrices:
X_temp = X;
Y_temp = Y;
t_matrix = ones(n, nlv);%Right loading vectors of X'*Y, weights for Y.
u_matrix = ones(n, nlv);%Left loading vectors of X'*Y, weights for X.
a_matrix = ones(p, nlv);
b_matrix = ones(q, nlv);
c_matrix = ones(p, nlv);%Loading vectors for variables of X.
d_matrix = ones(q, nlv);%Loading vectors for variables of Y in regression mode.
e_matrix = ones(q, nlv);%Loading vectors for variables of Y in canonical mode.
n_ones = ones(1, n);
p_ones = ones(1, p);
q_ones = ones(1, q);

%Missing values:
na_X = miss_values(X);
na_Y = miss_values(Y);
if sum(sum(na_X)) > 0
    isna_X = 1;
else
    isna_X = 0;
end
if sum(sum(na_Y)) > 0
    isna_Y = 1;
else
    isna_Y = 0;
end



%SPARSE PLS

%Loop for number of components:
for h = 1:nlv
    nx = p - nvarX(:, h);
    ny = q - nvarY(:, h);

    %Step a:
    X_1 = X_temp;
    Y_1 = Y_temp;

    %Replace missing values by zeros:
    if isna_X == 1 
        X_1(na_X == 1) = 0;
    end 
    if isna_Y==1 
        Y_1(na_Y == 1) = 0;
    end 
    
    %Step b: svd of crossproduct of X and Y based on simpls:
    M = X_1' * Y_1;%Crossproduct.
    if find(isnan(M)),
        break;
    end
    [U, ~, V] = svd(M, 0);
    a_old = U(:, 1);
    b_old = V(:, 1);
    if isna_X == 1 
        t = X_1 * a_old;
        A = repmat(a_old, size(n_ones));
        A(na_X' == 1) = 0;
        a_norm = A' * A;
        t = t ./ diag(a_norm);
        normt = sqrt(t' * t);
        t = t / normt;
    else
        t = X_temp * a_old / (a_old' * a_old);
        normt = sqrt(t' * t);
        t = t / normt;
    end
    if isna_Y == 1 
        u = Y_1 * b_old;
        B = repmat(b_old, size(n_ones));
        B(na_Y' == 1) = 0;
        b_norm = B' * B;
        u = u ./ diag(b_norm);
        normu = sqrt(u' * u);
        u = u / normu;
    else
        u = Y_temp * b_old / (b_old' * b_old);
        normu = sqrt(u' * u);
        u = u / normu;
    end   
    iter = 1;
    
    %Step c: loop for convergence of a and b:
    if isna_X == 1 
        a = X_1' * u;
    else 
        a = X_temp' * u;
    end
    if isna_Y == 1 
        b = Y_1' * t;
    else 
        b = Y_temp' * t;
    end
    d_a = a - a_old;
    d_b = b - b_old;
 
    while  d_a' * d_a > tol %|| d_b'*d_b > tol 
        if isna_X == 1 
            a = X_1' * u;
        else 
            a = X_temp' * u;
        end
        if isna_Y == 1 
            b = Y_1' * t;
        else 
            b = Y_temp' * t;
        end
        if nx ~= 0 
        	[a_sorted, ~]=sort(abs(a), 1);
            a_threshold = a_sorted(nx, 1);%The smallest coefficient to keep.
        	for z = 1 : size(a, 1)%Number of variables.
            	if abs(a(z, :)) > a_threshold
                	tmp = abs(a(z, :)) - a_threshold; %Modified GT 21/05/2014.
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
        	[b_sorted, ~] = sort(abs(b), 1);
            b_threshold = b_sorted(ny, 1);%The smallest coefficient to keep.
            for zz = 1:size(b, 1)%Number of variables.
              	if abs(b(zz, :)) > b_threshold
                  	tmp = abs(b(zz, :))- b_threshold;%Keep only variables with coefficients higher than threshold.
                   	b(zz, :) = tmp * sign(b(zz, :)); %Modified GT 21/05/2014.
                else
                    b(zz, :) = 0;
                end
            end
        end
        b = b / (t' * t);
        if isna_X == 1 
            t = X_1 * a;
            A = repmat(a, size(n_ones));
            A(na_X' == 1) = 0;
            a_norm = A' * A;
            t = t ./ diag(a_norm);
            normt = sqrt(t' * t);
            t = t / normt;
        else
            t = X_temp * a / (a' * a);
            normt = sqrt(t' * t);
            t = t / normt;
        end
        if isna_Y == 1 
            u = Y_1 * b;
            B = repmat(b, size(n_ones));
            B(na_Y' == 1) = 0;
            b_norm = B' * B;
            u = u ./ diag(b_norm);
            normu = sqrt(u' * u);
            u = u / normu;
        else
            u = Y_temp * b / (b' * b);
            normu = sqrt(u' * u);
            u = u / normu;
        end 

        iter = iter + 1;
        d_a = a - a_old;
        d_b = b - b_old;
        a_old = a;
        b_old = b;
        if iter > max_iter
            %disp('Maximum number of iteration is reached');
            break
        end 
    end
   
    %Step d: deflation of matrices:
    if isna_X == 1 
        X_1 = X_temp;
        X_1(na_X == 1) = 0;
        c = X_1' * t;
        T = repmat(t, size(p_ones));
        T(na_X' == 1) = 0;
        t_norm = T' * T;
        c = c ./ diag(t_norm);
    else
        c = X_temp' * t / (t' * t);
    end
    X_temp = X_temp - (t * c');
    if mode == 1
     if isna_Y == 1 %Regression mode.
          Y_1 = Y_temp;
        Y_1(na_Y == 1) = 0;
        d = Y_1' * t;
        T = repmat(Y_1, size(q_ones));
        T(na_Y' == 1) = 0;
        t_norm = T' * T;
        d = d ./ diag(t_norm);
     else
         d = Y_temp' * t / (t' * t);
     end
     Y_temp = Y_temp - (t * d');
     e = e_matrix(:, h);
    else 
        if isna_Y == 1 %Canonical mode.
        Y_1 = Y_temp;
        Y_1(na_Y == 1) = 0;
        e = Y_1' * u;
        U = repmat(Y_1, size(q_ones));
        U(na_Y' == 1) = 0;
        u_norm = U' * U;
        e = e ./ diag(u_norm);
        else
            e = Y_temp' * u / (u' * u);
         d = d_matrix(:, h);
        end
    end
    
     %Collection of results:
     t_matrix(:, h) = t;
     u_matrix(:, h) = u;
    
     a_matrix(:, h) = a;
     b_matrix(:, h) = b;
     c_matrix(:, h) = c;
     d_matrix(:, h) = d;
     e_matrix(:, h) = e;
end
 
%CALCULATE BETA
if mode == 1
    lY = d_matrix;
elseif mode == 2
    lY = e_matrix;
end
while isnan(rcond(((c_matrix(:, 1:nlv))' * a_matrix(:, 1:nlv))))
    nlv=nlv-1;
end
R = a_matrix(:, 1:nlv) / ((c_matrix(:, 1:nlv))' * a_matrix(:, 1:nlv));
A = (lY(:, 1:nlv))';

t_matrix(:, nlv+1:end) = [];
u_matrix(:, nlv+1:end) = [];

a_matrix(:, nlv+1:end) = [];
b_matrix(:, nlv+1:end) = [];
c_matrix(:, nlv+1:end) = [];
d_matrix(:, nlv+1:end) = [];
e_matrix(:, nlv+1:end) = [];

B = R * A;
B0 = mean_Y - mean_X * B;






%COLLECTION OF OUTPUT

%Old:
%model.c_matrix = c_matrix; %Loading vectors for variables of X.
%model.d_matrix = d_matrix; %Loading vectors for variables of Y in regression mode.
%model.e_matrix = e_matrix; %Loading vectors for variables of Y in canonical mode.
%model.input.keep_X = nvarX; %Input value for keep_X.
%model.input.keep_Y = nvarY; %Input value for keep_Y.
%model.input.ncomp = nlv; %Input value for nlv.
%model.input.X = X; %Input value for X.
%model.input.Y = Y; %Input value for Y.
%model.X_mean = mean_X; %Mean values for the variables in X.
%model.Y_mean = mean_Y; %Mean values for the variables in Y.
%model.X_scale = scale_X; %Standard deviations for the variables in X.
%model.Y_scale = scale_Y; %Standard deviations for the variables in Y.
%model.variates.X = t_matrix; %Including t_matrix: right loading vectors of X'*Y, weights for Y.
%model.variates.Y = u_matrix; %Including u_matrix: left loading vectors of X'*Y, weights for X.
%model.loadings.X = a_matrix; %Including a_matrix.
%model.loadings.Y = b_matrix; %Including b_matrix.
%model.input.mode = mode; %Input value for mode.



%New, edited by T. Offermans:

model.P = c_matrix; %X loadings.
if mode == 1
    model.Q = d_matrix; %Y loadings in regression mode.
elseif mode == 2
    model.Q = e.matrix; %Y loadings in canonical mode.
end
model.T = t_matrix; %X scores.
model.U = u_matrix; %Y Scores.
model.B = B; %Beta.
model.B0 = B0; %Intercept values.
model.R = R; %Weight vector.
model.input.X = X; %Input value for X.
model.input.Y = Y; %Input value for Y.
model.input.nlv = nlv; %Input value for nlv.
model.input.nvarX = nvarX; %Input value for nvarX.
model.input.nvarY = nvarY; %Input value for nvarY.
model.input.max_iter = max_iter; %Input value for max_iter.
model.input.tol = tol; %Input value for tol.
model.input.mode = mode; %Input value for mode.
model.input.mc = mc; %If mean-centering was performed or not.
model.mean.X = mean_X; %Mean values for the variables in X.
model.mean.Y = mean_Y; %Mean values for the variables in Y.
model.std.X = scale_X; %Standard deviations for the variables in X.
model.std.Y = scale_Y; %Standard deviations for the variables in Y.

%I do not fully understand what these are:
model.loadings.X = a_matrix; %Including a_matrix.
model.loadings.Y = b_matrix; %Including b_matrix.

model.q = size(Y, 2);
model.p = size(X, 2);

end




%SUBFUNCTIONS

function[nzv_ind, Xcorr] = nzv_elimination(X, n, p)
% Aim: eliminating near-zero variance predictors in X
% a predictor is classified as near-zero variance if the percentage of unique values in the samples is less than 10% 
threshold_nzv=round(n/10);
nzv=ones(1,p);
for i=1:p
    nzv(1,p)=length(unique(X(:,p)));
end
for i=1:p
    if nzv(1,p)<threshold_nzv
        nzv_ind(1,p)=1;
    else
        nzv_ind(1,p)=0;
    end
end
Xcorr=X(:,nzv_ind == 0);
end

%___________________________________
function [isna_X] =miss_values(data)
% Aim: finding indices of missing values
isna_X=zeros(size(data));
for i=1:size(data,1)
    isna_X(i,:) = isnan(data(i,:));
end
end
  