function table2latex(T,fileName)
  fid = fopen (fileName, "w");
  fprintf(fid,"\\begin{table} \n");
  numCol = size(T.var,2) + 1;
  fprintf(fid,"\\begin{tabular}{");
  fprintf(fid,repmat('l',[1 numCol]));
  fprintf(fid,"}\n");
  %fprintf("\hline \n");
  for ii = 1:size(T.mat,1)
    for jj = 1:size(T.mat,2)
      if jj == 1 && ii == 1
        fprintf(fid, " &");
        fprintf(fid, " ");
      elseif jj == 1
        fprintf(fid,strcat(T.source{ii-1,1}," &"));
        fprintf(fid," ");
      end

      if ii == 1 && jj == size(T.mat,2)
        fprintf(fid, T.var{1,jj});
        fprintf(fid," ");
      elseif ii == 1
        fprintf(fid, strcat(T.var{1,jj}," &"));
        fprintf(fid," ");
      elseif jj == size(T.mat,2)
        fprintf(fid, num2str(T.mat(ii,jj)));
        fprintf(fid," ");
      else
        fprintf(fid,strcat(num2str(T.mat(ii,jj))," &"));
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










