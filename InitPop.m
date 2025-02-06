function Chrom = InitPop(NIND, N)
%% 初始化种群函数
% 输入：
%   NIND - 种群大小，即个体的数量
%   N    - 每个个体的基因长度（对应于城市的数量）
% 输出：
%   Chrom - 初始化的种群矩阵，每一行表示一个个体的路径

Chrom = zeros(NIND, N);  % 预分配种群矩阵
for i = 1:NIND
    Chrom(i, :) = randperm(N);  % 为每个个体生成随机排列路径
end
end
