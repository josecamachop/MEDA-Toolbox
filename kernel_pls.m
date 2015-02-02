%1. F. Lindgren, P. Geladi and S. Wold, J. Chemometrics, 7, 45 (1993).
%4. S. De Jong and C. J. F. Ter Braak, J. Chemometrics, 8, 169 (1994).
%Dayal BS, MacGregor JF. Improved PLS algorithms. J. Chemom. 1997;11: 7385. 

% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 04/Jan/13.
%

function [beta,W,P,Q,R] = kernel_pls(XX,XY,A)
%XY=X'*Y; % compute the covariance
%XX=X'*X; % matrices
%A: pcs
M = size(XY,2);
W=[];
P=[];
Q=[];
R=[];
for i=1:A, % A=number of PLS components to be computed
    if M==1, % if there is a single response variable, compute the
        w=XY; % X-weights as shown here
    else % else
        [C,D]=eig(XY'*XY); % ?rst compute the eigenvectors of YTXXTX
        q=C(:,find(diag(D)==max(diag(D)))); %find the eigenvector corresponding to the largest eigenvalue
        w=(XY*q); % compute X-weights
    end
    w=w/sqrt(w'*w); % normalize w to unity
    r=w; % loop to compute ri
    for j=1:i-1,
        r=r-(P(:,j)'*w)*R(:,j);
    end
    tt=(r'*XX*r); % compute tTt
    p=(r'*XX)'/tt; % X-loadings
    q=(r'*XY)'/tt; % Y-loadings
    XY=XY-(p*q')*tt; % XTY de?ation
    W=[W w]; % storing loadings and weights
    P=[P p];
    Q=[Q q];
    R=[R r];
end
beta=R*Q';