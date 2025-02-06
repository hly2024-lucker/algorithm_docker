%% 计算每个个体的适应度值

% 输入：
% D 为路径矩阵的距离
% Chrom 为个体的染色体


function len = PathLength(D,Chrom)

[row,col]=size(D);

NIND = size(Chrom,1);

len = zeros(NIND,1);

for i=1:NIND

p=[Chrom(i,:) Chrom(i,1)];

i1=p(1:end-1);

i2=p(2:end);

len(i,1)=sum(D((i1-1)*col+i2)); %% 获取当前个体的路径长度，返回的是一个标量值


end
