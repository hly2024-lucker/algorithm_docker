function FitnV = Fitness(len)
%% 适应度函数
% 输入：
%   len - 各个体路径的长度向量 (TSP 问题中的距离)
% 输出：
%   FitnV - 各个体的适应度值

% 计算适应度 (路径越短，适应度越高)
FitnV = 1 ./ len;
end
