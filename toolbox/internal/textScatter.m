function ax = textScatter(fig_h,bdata,varargin)

% Print text in a Scatter plot.
%
% textScatter(fig_h,bdata) % minimum call
%
%
% INPUTS:
%
% fig_h: (1x1) figure handle
%
% bdata: (Nx2) bidimensional data 
%
%
% Optional INPUTS (parameter):
%
% 'EleLabel': [Nx1] name of the elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single
%   group by default)
%
% 'Multiplicity': [Nx1] multiplicity of each row (1s by default)
%
% 'PlotMult': str
%      'none': do not plot multiplicity (by default)
%      'zaxis': plot multiplicity information in the Z axis.
%      'zsize': plot multiplicity info in the size of the markers and
%               classes in Z-axis
%
% 'BlurIndex': [1x1] to avoid blur when adding labels. It reflects the
%   minimum distance (normalized to [0,1]) where a cluttered label is 
%   allowed to be visualized. For a value of 0, no cluttered labels are 
%   printed, while for a value of 1 all labels are printed, and thus the 
%   highest blur. By default 0.3 is chosen. 
%
%
% OUTPUTS:
%
% 'ax': [1x4] axis enclosing the text.
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 15/Jan/2025
%
% Copyright (C) 2025  University of Granada, Granada
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

figure(fig_h);

N = size(bdata, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',ones(N,1));  
addParameter(p,'PlotMult','none'); 
addParameter(p,'Multiplicity',ones(N,1)); 
addParameter(p,'BlurIndex',0.3);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
plottype = p.Results.PlotMult;
classes = p.Results.ObsClass;
mult = p.Results.Multiplicity;
blur = p.Results.BlurIndex;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel  = elabel';  end;
if size(classes,1) == 1, classes = classes'; end;
if size(mult,1) == 1, mult = mult'; end;

% Convert num arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end
if ~isempty(classes) && isnumeric(classes), classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;

% Validate dimensions of input data
assert(size(bdata,2) == 2, 'Dimension Error: paramter ''bdata'' must be N-by-2. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: paramter ''EleLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;


%% Main code

% Get ordering of classes
unique_classes = unique(classes,'stable');
if iscell(classes)
     ord_classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
     ord_classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end
unique_ord_classes = unique(ord_classes);

ax = axis;
deltax = (ax(2)-ax(1)); 
deltay = (ax(4)-ax(3)); 

c = 0.01; % bias with text
if ~isempty(elabel)
    for i=1:N
        suffx = length(char(strtrim(elabel(i,1))))+1;
        ind = [1:(i-1) (i+1):size(bdata,1)];
        
        dxM = (bdata(ind,1)-bdata(i,1))/(deltax*suffx/60); % app. 60 characters in the x-axis
        dxM(dxM<0) = Inf;
        dxM(dxM>1) = Inf;
        dyM = (bdata(ind,2)-bdata(i,2))/(deltay*2/40); % app. 40 characters in the y-axis
        dyM(dyM<-1) = Inf;
        dyM(dyM>1) = Inf;
        
        d = dxM.^2+dyM.^2;
        ind2 = find(min(d)==d,1);
        if dxM(ind2) > blur || dyM(ind2)==Inf  || isempty(ind)
            posx = bdata(i,1) + c*deltax;
            posy = bdata(i,2) + c*deltay;
            ax(2) = max(ax(2), posx + deltax*suffx/60);
            ax(4) = max(ax(4), posy + deltay*2/40);
            if strcmp(plottype,'zaxis')
                text(posx, posy, mult(i), strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 12);
            elseif strcmp(plottype,'zsize')
                text(posx, posy, ord_classes(i), strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 12);
            else
                text(posx, posy, strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 12);
            end
        else
            suffx(i) = 0;
        end
    end
end
