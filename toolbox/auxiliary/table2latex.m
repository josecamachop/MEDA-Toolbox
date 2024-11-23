function table2latex(T, fileName, formatspec, begend)

% table2latex converts the output of parglm.m (Parallel, General Linear
% Models) to a table formatted for use in LaTeX.
%
% table2latex(T, 'table.tex', '%.2f'); % minimum call
%
% See also: parglm
%
%
% INPUTS:
%
% T (table): MATLAB structure encoding a table as a numerical matrix with
% string headers.
%
% fileName (string): name of output file (e.g.: 'table.tex')
%
% formatspec (string): formatting of numerical values (e.g.: '%.2f')
%
% begend (int): toggles the header/footer characteristics of the output
% table, in case several tables are being grouped together as one. begend
% == 1: include \\begin{tabular}{ to initialize the tabular environment.
% begend == 2: include \\end{tabular} to terminate the tabular environment.
% begend == -1 include neither (for intermediate entries). begend == 0:
% include both \\begin{tabular} and \\end{tabular} for a standalone table.
%
% OUTPUTS:
%
% fileName.tex: a .tex file that can be input directly into LaTeX
%
%
% EXAMPLE OF USE (create table from parglm output)
%
% [parglmo,T] = parglm(X,F);
% table2latex(T,'outputTable.tex','%.2f',0)
%
% coded by: Michael Sorochan Armstrong (mdarmstr@ugr.es)
% last modification: 21/11/2024
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

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ~isOctave
    T2.source = table2cell(T(:,'Source'));
    T2.var = T.Properties.VariableNames;
    T2.mat = table2array(T(:,2:end));
    T = T2;
end

if nargin < 4
    begend = 0;
end

fid = fopen(fileName, 'w');

if begend == 0 || begend == 1
    fprintf(fid,'\\begin{tabular}{');
    fprintf(fid,repmat('l',[1 size(T.var,2) + 1]));
    fprintf(fid,'}\n');
end

if begend == 0 || begend == 1 || begend == 2 || begend == -1
    for ii = 1:size(T.mat,1)+1
        for jj = 1:size(T.mat,2)+1
            if jj == 1 && ii == 1
                fprintf(fid, ' &');
                fprintf(fid, ' ');
            elseif jj == 1
                fprintf(fid,strcat(T.source{ii-1,1},' &'));
                fprintf(fid,' ');
            elseif ii == 1 && jj == size(T.mat,2)+1
                fprintf(fid, T.var{1,jj});
                fprintf(fid,' ');
            elseif ii == 1
                fprintf(fid, strcat(T.var{1,jj},' &'));
                fprintf(fid,' ');
            elseif jj == size(T.mat,2)+1
                fprintf(fid, num2str(T.mat(ii-1,jj-1),formatspec));
                fprintf(fid,' ');
            else
                fprintf(fid,strcat(num2str(T.mat(ii-1,jj-1),formatspec),' &'));
                fprintf(fid,' ');
            end
        end

        if ii == 1
            fprintf(fid,'\\\\ \n');
            fprintf(fid,' \\hline \n');
        else
            fprintf(fid,'\\\\ \n');
        end
    end
end

if begend == 0 || begend == 2
    fprintf(fid,'\\end{tabular} \n');
end

fclose(fid);