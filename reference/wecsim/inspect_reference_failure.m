root = 'C:/Users/alex/topwater';
d = sllasterror;
file = fullfile(root,'artifacts/reference-imported-dispatch-stock/last-diagnostic.json');
fid = fopen(file,'w');
fprintf(fid,'%s\n',jsonencode(d,PrettyPrint=true));
fclose(fid);
fprintf('%s\n',jsonencode(d));
fprintf('Loaded models: %s\n',strjoin(find_system('type','block_diagram'),', '));
