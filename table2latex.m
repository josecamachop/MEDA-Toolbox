function table2latex(T,fileName)

% table2latex(T, "table.tex");
% Function to convert the table structure that is output from the MEDA toolbox to a latex markup for publication.
% Works with the Octave table structure, but assumes that the input is an Octave table structure if Octave is being used.
% Michael Sorochan Armstrong, 2022-12-19
%
% Create exception for edge case where Matlab is being used, but there is an Octave table structure.

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ~isOctave
  T2.source = table2cell(T(:,'Source'));
  T2.var = T.Properties.VariableNames;
  T2.mat = table2array(T(:,2:end))
  T = T2;
end

  fid = fopen (fileName, "w");
  fprintf(fid,"\\begin{table} \n");
  numCol = size(T.var,2) + 1;
  fprintf(fid,"\\begin{tabular}{");
  fprintf(fid,repmat('l',[1 numCol]));
  fprintf(fid,"}\n");
  for ii = 1:size(T.mat,1)+1
    for jj = 1:size(T.mat,2)+1
      if jj == 1 && ii == 1
        fprintf(fid, " &");
        fprintf(fid, " ");
      elseif jj == 1
        fprintf(fid,strcat(T.source{ii-1,1}," &"));
        fprintf(fid," ");
      elseif ii == 1 && jj == size(T.mat,2)+1
        fprintf(fid, T.var{1,jj-1});
        fprintf(fid," ");
      elseif ii == 1
        fprintf(fid, strcat(T.var{1,jj-1}," &"));
        fprintf(fid," ");
      elseif jj == size(T.mat,2)+1
        fprintf(fid, num2str(T.mat(ii-1,jj-1)));
        fprintf(fid," ");
      else
        fprintf(fid,strcat(num2str(T.mat(ii-1,jj-1))," &"));
        fprintf(fid," ");
      end
    end
    if ii == 1
      fprintf(fid,"\\\\ \n");
      fprintf(fid," \\hline \n");
    else
      fprintf(fid,"\\\\ \n");
    endif
  end
  fprintf(fid,"\\end{tabular} \n");
  fprintf(fid,"\\end{table} \n");
  fclose(fid);
end
