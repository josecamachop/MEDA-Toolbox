
function fig_h = plot_map(map,label)

% Plot color map.
%
% plot_map(map) % minimum call
% plot_map(map,label) % complete call
%
%
% INPUTS:
%
% map: (MxM) matrix with values in the [0,1] interval. 
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
%
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

if nargin < 1, error('Error in the number of arguments.'); end;
s = size(map);
if nargin < 2 || isempty(label)
    label=num2str((1:s(2))'); 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end

%% Main code

fig_h=figure;
map3 = [map map(:,end);map(end,:) map(end,end)];
sur_h=surface((1:s(2)+1)'*ones(1,s(2)+1),ones(s(2)+1,1)*(1:s(2)+1),map3);
axes_h = get(sur_h,'Parent');
set(axes_h,'Box','on');
set(axes_h,'XAxisLocation','top');
set(axes_h,'YDir','reverse');
set(axes_h,'XTick',(1:s(2))+0.5);
set(axes_h,'YTick',(1:s(2))+0.5);
set(axes_h,'XTickLabel',label);
set(axes_h,'YTickLabel',label);
set(axes_h,'FontSize',14);
ind = [0:.2:0.79 0.8:0.04:1]';
set(fig_h,'Colormap',[[ind;ones(10,1)] [ind;flipud(ind)] [ones(10,1);flipud(ind)]])
caxis([-1 1]);
if find(map>0 & map<1)
    c_h=colorbar;
    set(c_h,'FontSize',14);
end

axis([1,s(2)+1,1,s(2)+1]);

    


        