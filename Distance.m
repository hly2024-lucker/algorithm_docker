function D = Distance(a)
%% 计算城市之间的距离矩阵
% 输入: 
%   a - 每个城市的位置矩阵 (n x 2)，每行表示一个城市的 (x, y) 坐标
% 输出: 
%   D - 距离矩阵 (n x n)，D(i, j) 表示城市 i 到城市 j 的距离

n = size(a, 1);  % 获取城市数量
D = zeros(n, n); % 初始化距离矩阵

% 计算每两个城市之间的欧氏距离
for i = 1:n
    for j = 1:n
        D(i, j) = sqrt((a(i, 1) - a(j, 1))^2 + (a(i, 2) - a(j, 2))^2);
    end
end
end
