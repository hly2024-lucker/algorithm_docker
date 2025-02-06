clear
clc
close all

citys = load('lin318.txt');  % 加载城市数据
citys(:,1) = [];            % 去除第一列（假设不需要）
n = size(citys,1);          % 城市数量
D = zeros(n,n);             % 初始化距离矩阵
D = Distance(citys);        % 生成城市之间的距离矩阵

NIND = 45;  % 种群大小（蚂蚁数量）
alpha = 1;  % 信息素重要性 
beta = 5;   % 启发式因子重要性
rho = 0.15;  % 信息素挥发系数
Q = 1;      % 信息素常数
Eta = 1./D; % 启发式信息矩阵（距离的倒数）
Tau = ones(n,n);  % 初始化信息素矩阵
iter_max = 600;   % 最大迭代次数
N = size(D,1);          % 矩阵 D 的维度，即城市数量 (假设这里是 34x34)
m = NIND;               % 蚂蚁数量（种群大小）

Table = zeros(m,n);     % 路径记录表，每只蚂蚁的访问顺序（m 行 n 列）


Chrom = more_init(NIND, N, D);  % 初始化种群
Table = Chrom(1:NIND,:);        % 存储路径的表格
Length = PathLength(D, Chrom);  % 计算路径长度


iter = 1;                     % 迭代初始值
Route_best = zeros(iter_max, n);   % 每次迭代的最佳路径
Length_best = zeros(iter_max, 1);  % 每次迭代的最佳路径长度
Length_ave = zeros(iter_max, 1);   % 每次迭代的平均路径长度

%%%%%%%%%%%%%%%%%% 加入双向贪婪机制的初始化蚁群，并分配信息素 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 初始化种群信息素
% 逐个蚂蚁的路径计算
    Delta_Tau = zeros(n, n);  % 信息素增量矩阵
    k1=0;
    k2=0;
    min_len=min(Length);
    max_len=max(Length);
    
    %计算最差路径和最好路径的个数
    for i=1:m
        if min_len==Length(i)
            k1=k1+1;
        elseif max_len==Length(i)
            k2=k2+1;
        end
    end
    L_min=zeros(n,n);
    L_max=zeros(n,n);
    [min_Length,min_index]=min(Length);
    [max_Length,max_index]=max(Length);
    for i=1:(n-1)%局部最优最差路径
        L_min(Table(min_index,i),Table(min_index,i+1))=(Q/Length(min_index))*k1;
        L_max(Table(max_index,i),Table(max_index,i+1))=(Q/Length(max_index))*k2;
    end
    L_min(Table(min_index,n),Table(min_index,1))=(Q/Length(min_index))*k1;
    L_max(Table(max_index,n),Table(max_index,1))=(Q/Length(max_index))*k2;
    %(增强最好路径的信息素，忽略最差路径的信息素)
    Tau = (1 - rho) * Tau + Delta_Tau+L_min-L_max;  % 更新信息素矩阵，新的公式

%% 迭代寻找最优路径
while iter <= iter_max
    % 随机生成每只蚂蚁的起始城市
    start = zeros(m, 1);
    for i = 1:m
        temp = randperm(n);  % 随机排列所有城市编号
        start(i) = temp(1);  % 每只蚂蚁选择第一个城市作为起点
    end
    Table(:, 1) = start;  % 初始化路径记录表的第一列为起始城市

    % 构建空闲城市列表
    citys_index = 1:n;

    % 逐个蚂蚁选择路径
    for i = 1:m
        % 为每只蚂蚁选择完整路径
        for j = 2:n
            tabu = Table(i, 1:(j-1));  % 已访问的城市（禁忌列表）
            allow_index = ~ismember(citys_index, tabu);  % 剩余未访问城市
            allow = citys_index(allow_index);  % 可访问城市集合

            % 计算访问未访问城市的概率
            P = allow;
            for k = 1:length(allow)
                P(k) = Tau(tabu(end), allow(k))^alpha * Eta(tabu(end), allow(k))^beta;
            end
            P = P / sum(P);  % 归一化概率

            % 使用轮盘赌法选择下一个城市
            Pc = cumsum(P);  % 累加概率形成轮盘
            target_index = find(Pc >= rand, 1);  % 随机选择一个城市
            target = allow(target_index);  % 确定访问的城市
            Table(i, j) = target;  % 更新路径表
        end
    end

    % 计算每只蚂蚁的路径长度
    Length = zeros(m, 1);
    for i = 1:m
        Route = Table(i, :);
        for j = 1:(n - 1)
            Length(i) = Length(i) + D(Route(j), Route(j + 1));
        end
        Length(i) = Length(i) + D(Route(n), Route(1));  % 加上回到起点的距离
    end

    % 更新最优路径及其长度
    if iter == 1
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min_Length;
        Length_ave(iter) = mean(Length);
        Route_best(iter, :) = Table(min_index, :);
    else
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter - 1), min_Length);
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter, :) = Table(min_index, :);
        else
            Route_best(iter, :) = Route_best(iter - 1, :);
        end
    end

    % 更新信息素矩阵
    Delta_Tau = zeros(n, n);
 k1=0;
    k2=0;
    min_len=min(Length);
    max_len=max(Length);
    
    %计算最差路径和最好路径的个数
    for i=1:m
        if min_len==Length(i)
            k1=k1+1;
        elseif max_len==Length(i)
            k2=k2+1;
        end
    end
    L_min=zeros(n,n);
    L_max=zeros(n,n);
    [min_Length,min_index]=min(Length);
    [max_Length,max_index]=max(Length);
    for i=1:(n-1)%局部最优最差路径
        L_min(Table(min_index,i),Table(min_index,i+1))=(Q/Length(min_index))*k1;
        L_max(Table(max_index,i),Table(max_index,i+1))=(Q/Length(max_index))*k2;
    end
    L_min(Table(min_index,n),Table(min_index,1))=(Q/Length(min_index))*k1;
    L_max(Table(max_index,n),Table(max_index,1))=(Q/Length(max_index))*k2;
    %(增强最好路径的信息素，忽略最差路径的信息素)
    Tau = (1 - rho) * Tau + Delta_Tau+L_min-L_max;  % 更新信息素矩阵，新的公式


    % 迭代次数增加，清空路径表
    iter = iter + 1;
    Table = zeros(m, n);
end


%% 结果显示
[Shortest_Length, index] = min(Length_best);  % 找到最短路径及其索引
Shortest_Route = Route_best(index, :);  % 获取最优路径

% 打印最优路径及其长度
disp(['最短距离: ', num2str(Shortest_Length)]);
disp(['最优路径: ', num2str([Shortest_Route Shortest_Route(1)])]);  % 回到起点

%% 绘图
figure(1);
plot([citys(Shortest_Route, 1); citys(Shortest_Route(1), 1)], ...
     [citys(Shortest_Route, 2); citys(Shortest_Route(1), 2)], 'o-');  % 绘制路径
grid on;

% 为每个城市标注编号
for i = 1:size(citys, 1)
    text(citys(i, 1), citys(i, 2), ['  ' num2str(i)]);
end

% 标注起点和终点
text(citys(Shortest_Route(1), 1), citys(Shortest_Route(1), 2), ' 起点');
text(citys(Shortest_Route(end), 1), citys(Shortest_Route(end), 2), ' 终点');

% 设置坐标轴标签和标题
xlabel('城市位置 X 轴');
ylabel('城市位置 Y 轴');
title(['蚁群算法优化路径 (最短距离: ' num2str(Shortest_Length) ')']);

%% 路径长度变化曲线
figure(2);
plot(1:iter_max, Length_best, 'b', 1:iter_max, Length_ave, 'r');  % 绘制长度变化曲线
legend('最短距离', '平均距离');
xlabel('迭代次数');
ylabel('距离');
title('迭代过程中最短距离与平均距离变化');






