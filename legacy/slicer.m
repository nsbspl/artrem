%% parallelize: Take the data and cut it to 8 parts
n = 6; %number of threads to be generated
cutpoint = zeros (n+1,1);
for i = 2:n
    cutpoint(i) = ceil(i*dpts/n);
end
paralmatrix = zeros (ceil(dpts/n),2n);
for i = 1:n
    paral_matrix(2i-1:2i,i) = data_matrix(:,cutpoint(i):cutpoint(i+1));
end
1 
