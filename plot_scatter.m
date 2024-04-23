function fig_h = plot_scatter(bdata,varargin)

% Scatter plot.
%
% plot_scatter(bdata) % minimum call
% plot_scatter(bdata,elabel,classes,xylabel,lcont,0) % plot without multiplicity
% fig_h = plot_scatter(bdata,'EleLabel',elabel,'ObsClass',classes,'XYLabel',xylabel,'LimCont',lcont,'Option',opt,'Multiplicity',mult,'Threshold',maxv,'BlurIndex',blur) % complete call
%
%
% INPUTS:
%
% bdata: (Nx2) bidimensional data to plot.
%
% Optional INPUTS (parameters):
%
% 'EleLabel': [Nx1] name of the elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single
%   group by default)
%
% 'XYLabel': {2} xlabel and ylabel (nothing by default)
%
% 'LimCont': {2} control limits on x and y axis (nothing by default)
%
% 'Option': (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: plot for numerical classes (consistent with a colorbar)
%           1: plot for categorical classes (consistent with a legend)
%       b:
%           0: do not plot multiplicity
%           1: plot multiplicity
%       c: (for b 0)
%           0: filled marks
%           1: empty marks
%       c: (for b 1)
%           00: plot multiplicity info in the size of the markers.
%           01: plot multiplicity info in the form of the markers.
%           10: plot multiplicity information in the Z axis.
%           11: plot multiplicity info in the size of the markers and
%               classes in Z-axis
%
%   By deafult, opt = '100'. If less digits are specified, least significant
%   digits are set to 0, i.e. opt = 1 means a=1, b=0, c=0
%
% 'Multiplicity': [Nx1] multiplicity of each row (1s by default)
%
% 'Threshold': [1x3] thresholds for the different markers.
%       maxv(1): maximum threshold for marker 'd' for opt = 1101 (20 by default)
%       maxv(2): maximum threshold for marker 'o' for opt = 1101 (50 by default)
%       maxv(3): maximum threshold for marker 's' for opt = 1101 (100 by default)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels
%   are printer (the higher blur). Inf shows all the labels (1 by default).
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% EXAMPLE OF USE: Plot random data with filled marks and control limits:
%
% fig_h = plot_scatter(rand(100,2),'XYLabel',{'Y','X'},'LimCont',{0.8,0.8});
%
%
% EXAMPLE OF USE: with labels and classes in elements:
%
% fig_h = plot_scatter(randn(5,2),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'});
%
%
% EXAMPLE OF USE: with labels, multilicity and classes in elements:
%
% X = randn(5,2);
% opts = {'10' '1100' '1101' '1110' '1111'};
% for o = 1:length(opts),
%   plot_scatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'Option',opts{o},'Multiplicity',[1 20 50 100 1000]);
% end
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 23/Apr/2024
%
% Copyright (C) 2024  University of Granada, Granada
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

%

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

N = size(bdata, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',ones(N,1));
addParameter(p,'XYLabel',{'',''});
addParameter(p,'LimCont',[]);
addParameter(p,'Option','100');
addParameter(p,'Multiplicity',ones(N,1));
addParameter(p,'Threshold',[20 50 100]);
addParameter(p,'BlurIndex',1);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
classes = p.Results.ObsClass;
xylabel = p.Results.XYLabel;
lcont = p.Results.LimCont;
opt = p.Results.Option;
mult = p.Results.Multiplicity;
maxv = p.Results.Threshold;
blur = p.Results.BlurIndex;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel  = elabel';  end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;
if size(mult,1) == 1, mult = mult'; end;
if size(maxv,2) == 1, maxv = maxv'; end;

% Convert num arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Correct for opt integrity
if opt(1) == 0 && ~isnumeric(classes), opt(1) = 1; end
while length(opt)<4, opt = strcat(opt,'0'); end
if length(opt)<5 && opt(3)==1, opt = strcat(opt,'0'); end

% Convert num arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end
if ~isempty(classes) && isnumeric(classes) && opt(1)=='1', classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(xylabel),  xylabel = cellstr(xylabel); end;

% Validate dimensions of input data
assert(size(bdata,2) == 2, 'Dimension Error: parameter ''bdata'' must be N-by-2. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(elabel), assert((isequal(size(elabel), [N 1]) || isequal(size(elabel), [N+1 1])), 'Dimension Error: parameter ''EleLabel''  must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert ((isequal(size(classes), [N 1]) || isequal(size(classes), [N+1 1])), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: parameter ''XYLabel'' must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(lcont), assert (iscell(lcont) && isequal(size(lcont), [2 1]), 'Dimension Error: parameter ''LimCont'' must be a cell of 2 elements. Type ''help %s'' for more info.', routine(1).name); end;
assert (ischar(opt) && length(opt)==4, 'Dimension Error: parameter ''Option'' must be a string or num of maximum 4 bits. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(maxv), assert (isequal(size(maxv), [1 3]), 'Dimension Error: parameter ''Threshold'' must be 1-by-3. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Create figure window
fig_h = figure;
hold on;

% Sort data for colorbar
if opt(1)=='0'
    [classes,ord] = sort(classes,'ascend');
    cax = [min(classes) max(classes)];
    classes = num2str(classes); 
    classes = cellstr(classes);
    bdata = bdata(ord,:);
    elabel = elabel(ord);
    mult = mult(ord);
end

% Get ordering of classes
unique_classes = unique(classes,'stable');
if iscell(classes)
     ord_classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
     ord_classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end
unique_ord_classes = unique(ord_classes);

% Define mult bins, markers, colors and sizes
bins = [0 1 maxv Inf];
markers = ['^','v','d','o','s'];

if opt(1) == '0'
    color_list = parula(length(unique_ord_classes));
else
    color_list = hsv(length(unique_ord_classes));
end
    
colors = color_list(ord_classes, :);

sizes = zeros(size(mult));
for i=1:length(bins)-1
    sizes (mult>bins(i) & mult<=bins(i+1)) = round(2.5 * i^2 * pi);
end

switch opt(2:4)
    case '000'  % 2D plot, No multiplicity info, filled marks
        for i=1:length(unique_ord_classes)
            ind = find(ord_classes == unique_ord_classes(i));
            scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:),'filled','DisplayName',unique_classes{i});
        end

    case '010'  % 2D plot, No multiplicity info, empty marks
        for i=1:length(unique_ord_classes)
            ind = find(ord_classes == unique_ord_classes(i));
            scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:),'DisplayName',unique_classes{i});
        end

    case '100'  % 2D plot, Multiplicity in size
        for i=1:length(unique_ord_classes)
            ind = find(ord_classes == unique_ord_classes(i));
            scatter(bdata(ind,1), bdata(ind,2),sizes(ind), colors(ind,:),'filled','DisplayName',unique_classes{i});
        end

    case '101'  % 2D plot, Multiplicity in markers
        for i=1:length(unique_ord_classes)
            for j=1:length(bins)-1
                ind = find(ord_classes == unique_ord_classes(i) & mult<=bins(j+1) & mult>bins(j));
                disp_name = strcat(num2str(unique_ord_classes(i)), ' (mult: > ', num2str(bins(j)), ')');
                scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:), 'filled', markers(j), 'DisplayName', disp_name);
            end
        end

    case '110'  % 3D plot, Multiplicity in Z-axis
        for i=1:length(unique_ord_classes)
            ind = find(ord_classes == unique_ord_classes(i));
            scatter3(bdata(ind,1), bdata(ind,2), mult(ind), [], colors(ind,:), 'filled', 'DisplayName',unique_classes{i});
        end

    case '111'  % 3D plot, Multiplicity in size, classes in Z-axis
        for i=1:length(unique_ord_classes)
            ind = find(ord_classes == unique_ord_classes(i));
            scatter3(bdata(ind,1), bdata(ind,2), ord_classes(ind), sizes(ind), colors(ind,:), 'filled', 'DisplayName',unique_classes{i});
        end
end

text_scatter(fig_h,bdata,'EleLabel',elabel,'ObsClass',classes,'Option',opt(2:4),'Multiplicity',mult,'BlurIndex',blur);

ax = axis;
ax([1 3]) = min(ax([1 3]),zeros(1,2));
ax([2 4]) = max(ax([2 4]),zeros(1,2));

if ~isempty(lcont) % Plot control limits
    if ~isempty(lcont{2})
        ax(3) = min([ax(3);lcont{2}(:)]);
        ax(4) = max([ax(4);lcont{2}(:)]);
    end
    if ~isempty(lcont{1})
        ax(1) = min([ax(1);lcont{1}(:)]);
        ax(2) = max([ax(2);lcont{1}(:)]);
        for i=1:length(lcont{1})
            plot([lcont{1}(i) lcont{1}(i)], ax(3:4), 'r--','LineWidth',2, 'HandleVisibility', 'off');
        end
    end
    if ~isempty(lcont{2})
        for i=1:length(lcont{2})
            plot(ax(1:2),[lcont{2}(i) lcont{2}(i)], 'r--','LineWidth',2, 'HandleVisibility', 'off');
        end
    end
else % Plot origin lines
    plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
    plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');
end
axis(ax)

% Set axis labels
if ~isempty(xylabel)
    xlabel(xylabel{1}, 'FontSize', 20);
    ylabel(xylabel{2}, 'FontSize', 20);
end

axes_h = get(fig_h,'Children');
for i=1:length(axes_h)
    if strcmp(get(axes_h(i), 'type'), 'axes')
        set(axes_h(i), 'FontSize', 12);
    end
end

% Set caxis if colorbar
if opt(1)=='0'
    caxis(cax);
end


legend off
box on
hold off



