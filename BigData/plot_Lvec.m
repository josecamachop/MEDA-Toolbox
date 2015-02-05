
function fig_h = plot_Lvec(vec,mult,olabel,slabel,lcont,maxv)

% Stem plot for large data.
%
% plot_vec(vec) % minimum call
% plot_Lvec(vec,olabel,slabel,lcont,maxv) % complete call
%
%
% INPUTS:
%
% vec: (Mx1) vector to plot. 
%
% mult: (Nx1) multiplicity of each observation (1s by default)
%
% olabel: (Mx1) name of the x-variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
% slabel: (str) y-variable/statistic plotted (nothing by default)
%
% lcont: (2xM) control limits.
%
% maxv: (1x3) thresholds for the different markers.
%       maxv(1): maximum threshold for marker 'x' for opt = 2 (20 by default)
%       maxv(1): maximum threshold for marker 'o' for opt = 2 (50 by default)
%       maxv(1): maximum threshold for marker 's' for opt = 2 (100 by default)
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 07/May/13.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

if nargin < 1, error('Error in the number of arguments.'); end;
s = max(size(vec));
if nargin < 2, mult = ones(s,1); end;
if nargin < 6, maxv = [20 50 100]; end;
if ndims(maxv)==2 & find(size(maxv)==max(size(maxv)))==1, maxv = maxv'; end;
if size(maxv,2)~=3 || size(maxv,1)~=1, error('Error in the dimension of the arguments.'); end;

%% Main code

maxv = [0 1 maxv Inf];

fig_h=figure;
hold on;
for j=1:length(maxv)-1,
    ind2 = find(mult<=maxv(j+1));
    ind3 = find(mult(ind2)>maxv(j));
    ind2 = ind2(ind3);
    stem(ind2,vec(ind2),'MarkerSize',2*j);
end

axes_h=get(fig_h,'Children');
axes_h=axes_h(1);
if exist('olabel')&~isempty(olabel), 
    set(axes_h,'XTick',1:s);
    set(axes_h,'XTickLabel',olabel);
end
if exist('slabel')
    ylabel(slabel,'FontSize',16);
end
set(axes_h,'FontSize',14);
if exist('lcont') && ~isempty(lcont)
    hold on
    b = [0.5:(s+1);0.5:(s+1)];
    a = [lcont(1,:);lcont(1,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2);
    a = [lcont(2,:);lcont(2,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2);
end    
axis tight
ax=axis;
axis auto
ax2=axis;
axis([ax(1:2) ax2(3:4)])
    


        