function DrawPath(Chrom, X)
%% 绘制路径的函数
% 输入:
%   Chrom - 表示路径的矩阵，每一行是一条路径
%   X     - 每个城市的坐标位置矩阵 (n x 2)，每行表示一个城市的 (x, y) 坐标

R = [Chrom(1, :) Chrom(1, 1)];  % 封闭路径，将起点添加到末尾
figure;
hold on;

% 绘制所有城市的坐标点
plot(X(:, 1), X(:, 2), 'o', 'Color', [0.5, 0.5, 0.5]);  
% 标记起点城市
plot(X(Chrom(1, 1), 1), X(Chrom(1, 1), 2), 'rv', 'MarkerSize', 20); 

% 标注每个城市的编号
for i = 1:size(X, 1)
    text(X(i, 1) + 0.05, X(i, 2) + 0.05, num2str(i), 'Color', [1, 0, 0]);
end

% 绘制路径并添加箭头
A = X(R, :);  % 提取路径上的城市坐标
for i = 2:size(A, 1)
    % 坐标转换到 figure 坐标系
    [arrowx, arrowy] = dsxy2figxy(gca, A(i-1:i, 1), A(i-1:i, 2));
    % 绘制带箭头的路径
    annotation('textarrow', arrowx, arrowy, 'HeadWidth', 8, 'Color', [0, 0, 1]);
end

hold off;
xlabel('横坐标');
ylabel('纵坐标');
title('路径图');
box on;

end
