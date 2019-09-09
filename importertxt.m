%% ImporterTxt
fid = fopen('chenlab_file1_data.txt');
data = textscan(fid, '%*s %f %*[^\n]','HeaderLines',1);
fid = fclose(fid);
vec = data{1,1};
L = length(vec);
data_Matrix = zeros(L,NumCol);
for i = 1:NumCol
fid = fopen('chenlab_file1_data.txt');
data = textscan(fid, [repmat('%*s',1,i-1), '%f', '%*[^\n]'],'HeaderLines',1);
fid = fclose(fid); 
data_Matrix(:,i) = data{1,1};
end