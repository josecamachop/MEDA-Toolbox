
function fig_h = plot_map(map,label,int,ind)

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
% label: (Mx1) name of the variables (numbers are used by default)
%
% int: (2x1) color interval ([-1;1] by default)
%
% ind: (Lx1) color distribution ([0:.2:0.79 0.8:0.04:1]' by default);
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 06/Apr/16.
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
M = size(map,2);
if nargin < 2 || isempty(label), label= []; end
if nargin < 3 || isempty(int), int = [-1;1]; end;
if nargin < 4 || isempty(ind), ind = [0:.2:0.79 0.8:0.04:1]'; end;

% Convert row arrays to column arrays
if size(label,1)  == 1, label = label'; end;

% Convert int arrays to str
if ~isempty(label) && isnumeric(label),
    vecn = label; 
    veci = 1:length(vecn);
    if length(label)>2, 
        max_lab = 30; % limit the number of labels displayed
        ini = 2;
        stepN = [];
        while isempty(stepN),
            lenv = length(vecn(ini:end));
            div = 1:(lenv-1);
            div = div(rem(lenv,div)==0);
            stepN = div(find(div>lenv/max_lab,1));
            ini = ini+1;
        end
        veci = 1:(lenv+ini-2);
        veci = veci(round([1 (ini-2+stepN):stepN:end]));
    end
    for i=veci,
        labele{i} = num2str(vecn(i));
    end
    label=labele'; 
end

% Convert char arrays to cell
if ischar(label),  label = cellstr(label); end;

%% Main code
fig_h=figure;
map3 = [map map(:,end);map(end,:) map(end,end)];
sur_h=surface((1:M+1)'*ones(1,M+1),ones(M+1,1)*(1:M+1),map3);
if M < 100
    set(sur_h,'EdgeColor',[0.95 0.95 0.95]);
else
    set(sur_h,'EdgeColor','none');
end

% Label font size
axes_h = get(sur_h,'Parent');
if ~isempty(label)
    label_length = max(cellfun('length', label));
    label_size = 300/(length(find(~cellfun('isempty', label)))*label_length);
    set(axes_h, 'FontSize', max(min(14,round(label_size)), 10));
end

% Set axis properties
set(axes_h,'Box','on');
set(axes_h,'XAxisLocation','top');
set(axes_h,'YDir','reverse');
if ~isempty(label)
    stepY = ceil(0.05*M/label_size);
    stepX = ceil(0.2*M/label_size);
    if stepX==1,
        set(axes_h,'XTick',(1:M)+0.5);
        set(axes_h,'XTickLabel',label);
    else
        set(axes_h,'XTickLabel','');
    end
    if stepY==1,
        set(axes_h,'YTick',(1:M)+0.5);
        set(axes_h,'YTickLabel',label);
    else
        set(axes_h,'YTickMode','auto');
    end
end

% Resize axes position
pos = get(axes_h, 'Position');
set(axes_h,'Position',[pos(1) pos(2)/2 pos(3) pos(4)])

% Set colors
set(fig_h,'Colormap',[[ind;ones(length(ind),1)] [ind;flipud(ind)] [ones(length(ind),1);flipud(ind)]])
caxis(int);
if find(map>0 && map<1)
    c_h=colorbar;
    set(c_h,'FontSize',14);
end

axis([1,M+1,1,M+1]);

    


        