function [Results_sparsePLS]=sparsepls(X,Y,ncomp, keepX, keepY, max_iter, tol, mode)
% Function to calculate sparse PLS models based on spls function by
% Sebastien Dejuan, France written in R, it uses procedure described in Le
% Cao et al.: A Sparse PLS for variable selection when intergrating Omics
% data, 2009
% Function handels missing values
% input: X, Y matrices, 
%        ncomp - number of latent variables, 
%        max_iter - maximum number of iteration e.g. 500, 
%        tol, - convergence criteria e.g. 1e-06
%        keepX, keepY - column vector indicating number of variables kept
%        per latent variable
%        mode: 1 - regression, 2 - canonical
% output: Results_sparsePLS:
    % c_matrix;loading vectors for variables of X
    % d_matrix; loading vectors for variables of Y in regression mode
    % e_matrix;loading vectors for variables of Y in canonical mode
    % input.: 
        % keepX;
        % keepY;
        % ncomp;
        % X;
        % Y;
        % mode: 1 - regression mode, 2 - canonical mode
    % mean_X;
    % mean_Y;
    % scale_X;
    % scale_Y;
    % variates.X including t_matrix: right loading vectors of X'*Y, weights for Y
    % variates.Y including u_matrix;  left loading vectors of X'*Y, weights for X
    % loadings.X including a_matrix;
    % loadings.Y including b_matrix;
    
% 01-12/05/2014 written and adopted by E.Szymanska, Radboud University Nijmegen
% 21/05/2014 edited by G.H. Tinnevelt, Radboud University Nijmegen
% VALIDATION of INPUT
% check size of X and Y matrix
[n,p]=size(X);
[r,q]=size(Y);
if length(size(X))~=2 % works only for matrices with more dimensions
    disp ('X must be a numeric matrix');
    return
end
if n ~=r
     disp ('Unequal number of rows in X and Y');
     return
end
% check other input parameters
if ncomp <= 0
    disp ('Invalid number of latent variables:');
    disp (ncomp);
    return
end
if ncomp > p
    warning('Maximum number of latent variables is set');
    ncomp=p;
end
if length(keepX) ~= ncomp
    disp ('Length of keepX must be equal to');
    disp (ncomp);
    return
end
if length(keepY) ~= ncomp
    disp ('Length of keepY must be equal to');
    disp (ncomp);
    return
end
if any(keepX > p)
    disp ('Each component of keepX must be lower or equal than');
    disp (p);
    return
end
if any(keepY > q)
    disp ('Each component of keepY must be lower or equal than');
    disp (q);
    return
end   
% setting default input parameters if they are not given
if nargin <6 
    max_iter=500;
    tol=0.000001;
    mode = 1;
    %keepX=repmat(100,1,ncomp);
    %keepY=repmat(100,1,ncomp);
end
% eliminating near-zero variance predictors in X
[nzv_ind, X]=nzv_elimination(X,n,p);    
% INITIALIZATION
% variable names
var_names_X=1:p;
var_names_Y=1:q;
sample_names=1:n;
% centring
mean_Y=mean(Y);
mean_X=mean(X);
scale_Y=(std(Y));
scale_X=(std(X));
Y=Y-repmat(mean(Y),n,1);
X=X-repmat(mean(X),n,1);
% autoscaling (in original code)
Y=Y./repmat(scale_Y, n, 1); 
X=X./repmat(scale_X, n, 1);
% creating intial matrices
X_temp=X;
Y_temp=Y;
t_matrix=ones(n,ncomp);% right loading vectors of X'*Y, weights for Y
u_matrix=ones(n,ncomp);% left loading vectors of X'*Y, weights for X
a_matrix=ones(p,ncomp);
b_matrix=ones(q,ncomp);
c_matrix=ones(p,ncomp);% loading vectors for variables of X
d_matrix=ones(q,ncomp);% loading vectors for variables of Y in regression mode
e_matrix=ones(q,ncomp);% loading vectors for variables of Y in canonical mode
n_ones=ones(1,n);
p_ones=ones(1,p);
q_ones=ones(1,q);
% missing values
na_X = miss_values(X);
na_Y = miss_values(Y);
if sum(sum(na_X)) >0
    isna_X=1;
else isna_X=0;
end
if sum(sum(na_Y)) >0
    isna_Y=1;
else isna_Y=0;  
end
% SPARSE PLS
% loop for number of components
for h=1:ncomp
    nx=p-keepX(:,h);
    ny=q-keepY(:,h);
    % step a 
    X_1=X_temp;
    Y_1=Y_temp;
    % missing values replacement by 0's
    if isna_X==1 
        X_1(na_X==1)=0;
    end 
    if isna_Y==1 
        Y_1(na_Y==1)=0;
    end 
    % step b: svd of crossproduct of X and Y based on simpls
    M=X_1'*Y_1;%crossproduct
    [U,~,V]=svd(M);
    a_old=U(:,1);
    b_old=V(:,1);
    if isna_X==1 
        t=X_1*a_old;
        A=repmat(a_old,size(n_ones));
        A(na_X'==1)=0;
        a_norm=A'*A;
        t=t./diag(a_norm);
        normt=sqrt(t'*t);
        t=t/normt;
    else t=X_temp*a_old/(a_old'*a_old);
        normt=sqrt(t'*t);
        t=t/normt;
    end
    if isna_Y==1 
        u=Y_1*b_old;
        B=repmat(b_old,size(n_ones));
        B(na_Y'==1)=0;
        b_norm=B'*B;
        u=u./diag(b_norm);
        normu=sqrt(u'*u);
        u=u/normu;
    else u=Y_temp*b_old/(b_old'*b_old);
        normu=sqrt(u'*u);
        u=u/normu;
    end   
    iter=1;
    % step c: loop for convergence of a and b
    if isna_X==1 
        a = X_1'*u;
    else 
        a= X_temp'*u;
    end
    if isna_Y==1 
        b = Y_1'*t;
    else 
        b= Y_temp'*t;
    end
    d_a=a-a_old;
    d_b=b-b_old;
 
    while  d_a'*d_a > tol %|| d_b'*d_b > tol 
        if isna_X==1 
            a = X_1'*u;
        else 
            a= X_temp'*u;
        end
        if isna_Y==1 
            b = Y_1'*t;
        else 
            b= Y_temp'*t;
        end
        if nx ~=0 
            [a_sorted,~]=sort(abs(a),1);
            a_threshold=a_sorted(nx,1);% the smallest coefficient to keep
           for z=1:size(a,1)% number of variables
               if abs(a(z,:)) > a_threshold
                   tmp=abs(a(z,:)) - a_threshold; %modified GT 21/05/2014
                   a(z,:) = tmp*sign(a(z,:));% keep only variables with coefficients higher than threshold
               else a(z,:)=0;
               end
           end
        end
        %a = a/(u'*u);
       norma=sqrt(a'*a);
       a=a/norma;
         if ny ~=0
             [b_sorted,~]=sort(abs(b),1);
                b_threshold=b_sorted(ny,1);% the smallest coefficient to keep
               for zz=1:size(b,1)% number of variables
                   if abs(b(zz,:))> b_threshold
                       tmp = abs(b(zz,:))- b_threshold;% keep only variables with coefficients higher than threshold
                       b(zz,:) = tmp*sign(b(zz,:)); %modified GT 21/05/2014
                   else b(zz,:)=0;
                   end
               end
         end
           b=b/(t'*t);
        if isna_X==1 
            t=X_1*a;
            A=repmat(a,size(n_ones));
            A(na_X'==1)=0;
            a_norm=A'*A;
            t=t./diag(a_norm);
            normt=sqrt(t'*t);
            t=t/normt;
        else
            t=X_temp*a/(a'*a);
            normt=sqrt(t'*t);
            t=t/normt;
        end
        if isna_Y==1 
            u=Y_1*b;
            B=repmat(b,size(n_ones));
            B(na_Y'==1)=0;
            b_norm=B'*B;
            u=u./diag(b_norm);
            normu=sqrt(u'*u);
            u=u/normu;
        else
            u=Y_temp*b/(b'*b);
            normu=sqrt(u'*u);
            u=u/normu;
        end 

        iter=iter+1;
        d_a=a-a_old;
        d_b=b-b_old;
        a_old=a;
        b_old=b;
        if iter > max_iter
            disp('Maximum number of iteration is reached');
            break
        end 
    end
   
    % step d: deflation of matrices
    if isna_X==1 
        X_1=X_temp;
        X_1(na_X==1)=0;
        c=X_1'*t;
        T=repmat(t,size(p_ones));
        T(na_X'==1)=0;
        t_norm=T'*T;
        c=c./diag(t_norm);
    else c=X_temp'*t/(t'*t);
    end
    X_temp=X_temp-(t*c');
    if mode==1
     if isna_Y==1 % regression mode
          Y_1=Y_temp;
        Y_1(na_Y==1)=0;
        d=Y_1'*t;
        T=repmat(Y_1,size(q_ones));
        T(na_Y'==1)=0;
        t_norm=T'*T;
        d=d./diag(t_norm);
     else d=Y_temp'*t/(t'*t);
     end
     Y_temp=Y_temp-(t*d');
     e=e_matrix(:,h);
    else 
        if isna_Y==1 % canonical mode 
          Y_1=Y_temp;
        Y_1(na_Y==1)=0;
        e=Y_1'*u;
        U=repmat(Y_1,size(q_ones));
        U(na_Y'==1)=0;
        u_norm=U'*U;
        e=e./diag(u_norm);
     else e=Y_temp'*u/(u'*u);
         d=d_matrix(:,h);
        end
    end
     % collection of results
     t_matrix(:,h)=t;
     u_matrix(:,h)=u;
     a_matrix(:,h)=a;
     b_matrix(:,h)=b;
     c_matrix(:,h)=c;
     d_matrix(:,h)=d;
     e_matrix(:,h)=e;
end % for loop for h
        
% COLLECTION OF OUTPUT
Results_sparsePLS.c_matrix=c_matrix;
Results_sparsePLS.d_matrix=d_matrix;
Results_sparsePLS.e_matrix=e_matrix;
Results_sparsePLS.input.keep_X=keepX;
Results_sparsePLS.input.keep_Y=keepY;
Results_sparsePLS.input.ncomp=ncomp;
Results_sparsePLS.input.X=X;
Results_sparsePLS.input.Y=Y;
Results_sparsePLS.X_mean=mean_X;
Results_sparsePLS.Y_mean=mean_Y;
Results_sparsePLS.X_scale=scale_X;
Results_sparsePLS.Y_scale=scale_Y;
Results_sparsePLS.variates.X=t_matrix;
Results_sparsePLS.variates.Y=u_matrix;
Results_sparsePLS.loadings.X=a_matrix;
Results_sparsePLS.loadings.Y=b_matrix;
Results_sparsePLS.input.mode=mode;
end
% SUBFUNCTIONS
%__________________________________________ 

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
    else nzv_ind(1,p)=0; 
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
  