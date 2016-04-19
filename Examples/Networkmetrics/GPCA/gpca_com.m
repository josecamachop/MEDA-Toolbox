% EDA example with GPCA using the MEDA Toolbox. See README.txt for more
% details.
% 
%  Data set and Analysis: 
% 
%    J. Camacho, R.A. Rodríguez-Gómez, E. Saccenti, Group-wise Principal 
% Component Analysis for Exploratory Data Analysis. Submitted to JCGS.

% coded by: José Camacho Páez.
% last modification: 19/Apr/16.

%% Inicialization, remember to set the path of the toolbox

clear
close all
clc

load gpca

%% Selection of the PCs

var_pca(x,0:30,2); % 15 PCs could be a choice

%% Visualize MEDA

[meda_map,meda_dis,ord] = meda_pca(x,1:15,2,0.1,11); % several groups of variables are found


%% Compute states with GIA

c = 0.7;

[b ,stg]=  gia(meda_map,c);


%% Prepare data for treemap visualization: load fortreemap to http://nesg.ugr.es/meda-visualization/

names = var_l;

weights = (weight-1)/9;

states = b';
    
save fortreemap names weights states


%% G-PCA

xcs = preprocess2D(x,2);

[P,T] =  gpca(xcs,stg,1:2);

%% Visualizing scores and loads

for i=1:2,
    t = T(:,i);
    lim = tinv(1-0.99,size(t,1)-1)*std(t)*sqrt(1+1/size(t,1));
    figure,
    subplot(2,1,1), plot(t,'.-'), hold on, plot(-lim*ones(size(t)),'r--'),plot(lim*ones(size(t)),'r--')
    subplot(2,1,2), bar(P(ord,i))
end


