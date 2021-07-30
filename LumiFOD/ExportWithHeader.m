function fileID=ExportWithHeader(fname,header,num_data,track_info)
% generate a new file and pointer
fileID = fopen(fname,'w');
% write the header as comma separated 
formatSpec=[];
for i=1:numel(header)
    tmp=header{i};
    fprintf(fileID,' %s ,', tmp);
end
fprintf(fileID,'\n');
% save the numerical data
for i=1:size(num_data,1)
    % write track_info
    fprintf(fileID,' %s ,',track_info{i});
    % write stats
    for j=1:size(num_data,2)
        fprintf(fileID,'%0.4f ,',num_data(i,j));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);