
function fig_h = plot_map(map,varargin)

% Plot color map.
%
% plot_map(map) % minimum call
% plot_map(map,'VarsLabel',label,'ColorInt',int,'ColorDist',ind) % complete call
%
%
% INPUTS:
%
% map: (MxM) matrix with values in the [0,1] interval. 
%
% Optional INPUTS (parameters):
%
% 'VarsLabel': (Mx1) name of the variables (numbers are used by default)
%
% 'ColorInt': (2x1) color interval ([-1;1] by default)
%
% 'ColoDist': (Lx1) color distribution ([0:.2:0.79 0.8:0.04:1]' by default);
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% plot_map(corr(X));
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 23/Apr/2024
%
% Copyright (C) 2024  University of Granada, Granada
% Copyright (C) 2024  Jose Camacho Paez, Alejandro Perez Villegas
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

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'VarsLabel',1:M);   
addParameter(p,'ColorInt',[-1;1]);
addParameter(p,'ColorDist',[0:.2:0.79 0.8:0.04:1]');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
label = p.Results.VarsLabel;
int = p.Results.ColorInt;
ind = p.Results.ColorDist;


% Convert row arrays to column arrays
if size(label,1)  == 1, label = label'; end;

% Convert int arrays to str
if ~isempty(label) && isnumeric(label)
    vecn = label; 
    veci = 1:length(vecn);
    if length(label)>2 
        max_lab = 30; % limit the number of labels displayed
        ini = 2;
        stepN = [];
        while isempty(stepN)
            lenv = length(vecn(ini:end));
            div = 1:(lenv-1);
            div = div(rem(lenv,div)==0);
            stepN = div(find(div>lenv/max_lab,1));
            ini = ini+1;
        end
        veci = 1:(lenv+ini-2);
        veci = veci(round([1 (ini-2+stepN):stepN:end]));
    end
    for i=veci
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
    lablength = cellfun('length', label);
    label_length = max(lablength(1:end-1)+lablength(2:end))/2;
    label_sizeH = 5/label_length;
    label_sizeV = 5;
    set(axes_h, 'FontSize', max(min(14,round(label_sizeH)), 10));
end

% Set axis properties
set(axes_h,'Box','on');
set(axes_h,'XAxisLocation','top');
set(axes_h,'YDir','reverse');
if ~isempty(label) 
    
    MaxRot = 60;
    if 0.05*M<label_sizeH % labels do not need to be rotated
        vals = 1:M;
        set(axes_h,'XTick',vals + 0.5);
        set(axes_h,'XTickLabel',label(vals));
    elseif 0.05*M<label_sizeV % labels are rotated
        vals = 1:M;
        set(axes_h,'XTick',vals + 0.5);
        set(axes_h,'XTickLabel',label(vals));
        set(axes_h,'XTickLabelRotation',ceil(MaxRot*0.05*M/label_sizeV));
    else % labels are reduced
        set(axes_h,'XTickMode','auto');
        ind1 = get(axes_h,'XTick');
        ind2 = find(ind1>0&ind1<=length(label));
        set(axes_h,'XTick',ind1(ind2) + 0.5);
        set(axes_h,'XTickLabel',label(ind1(ind2)));
        set(axes_h,'XTickLabelRotation',MaxRot);
        set(axes_h, 'FontSize', 14);
    end
    
    if 0.05*M<label_sizeV % labels do not need to be rotated
        vals = 1:M;
        set(axes_h,'YTick',vals + 0.5);
        set(axes_h,'YTickLabel',label(vals));
    else % labels are reduced
        set(axes_h,'YTickMode','auto');
        ind1 = get(axes_h,'YTick');
        ind2 = find(ind1>0&ind1<=length(label));
        set(axes_h,'YTick',ind1(ind2) + 0.5);
        set(axes_h,'YTickLabel',label(ind1(ind2)));
        set(axes_h,'FontSize', 14);
    end

end

% Resize axes position
pos = get(axes_h, 'Position');
set(axes_h,'Position',[pos(1) pos(2)/2 pos(3) pos(4)])

% Set colors
if int(1)<0
    set(fig_h,'Colormap',[[ind;ones(length(ind),1)] [ind;flipud(ind)] [ones(length(ind),1);flipud(ind)]])
else
    set(fig_h,'Colormap',[[ones(length(ind),1)] [flipud(ind)] [flipud(ind)]])
end
caxis(int);
if find(map>0 & map<=1)
    c_h=colorbar;
    set(c_h,'FontSize',14);
end

axis([1,M+1,1,M+1]);

    


        