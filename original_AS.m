clear
clc
close all

citys = load('bier127.txt');  % 加载城市坐标数据
citys(:,1) = [];            % 移除第一列（假设是城市编号）
n = size(citys,1);          % 城市数量
D = zeros(n,n);             
NIND = 45;                  % 种群大小
D = Distance(citys);        % 生成距离矩阵
N = size(D,1);              % 获取城市数量（假设为34）

m = 45;                     % 蚂蚁数量
alpha = 1;                  % 信息素重要性因子
beta = 5;                   % 启发因子重要性权重
rho = 0.15;                  % 信息素挥发因子
Q = 1;                      % 常系数
Eta = 1./D;                 % 启发函数
Tau = ones(n,n);            % 信息素矩阵
Table = zeros(m,n);         % 路径记录表
iter = 1;                   % 迭代初始值
iter_max = 400;             % 最大迭代次数
Route_best = zeros(iter_max,n);     % 每代最优路径
Length_best = zeros(iter_max,1);    % 每代最优路径长度
Length_ave = zeros(iter_max,1);     % 每代路径的平均长度


%% 迭代查找最优路径
while iter <= iter_max
    % 随机生成每个蚂蚁的起点城市
    start = zeros(m,1);
    for i = 1:m
        temp = randperm(n);
        start(i) = temp(1);
    end
    Table(:,1) = start;
    
    % 创建一个空白空间
    citys_index = 1:n;
    
    % 每只蚂蚁路径选择
    for i = 1:m
        % 每个城市路径选择
        for j = 2:n
            tabu = Table(i,1:(j-1));    % 已访问的城市集合（禁忌表）
            allow_index = ~ismember(citys_index,tabu);
            allow = citys_index(allow_index);   % 可访问的城市集合
            P = allow;
            
            % 计算城市间转移概率
            for k = 1:length(allow)
                P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
            end
            P = P / sum(P);
            
            % 轮盘赌法选择下一个访问城市
            Pc = cumsum(P); % 累加
            target_index = find(Pc >= rand);
            target = allow(target_index(1));
            Table(i,j) = target;
        end
    end

    Length = zeros(m,1);
    for i = 1:m
        Route = Table(i,:);
        for j = 1:(n-1)
            Length(i) = Length(i) + D(Route(j), Route(j+1));
        end
        Length(i) = Length(i) + D(Route(n), Route(1));
    end

    % 计算最短路径距离及平均距离
    if iter == 1
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min_Length;
        Length_ave(iter) = mean(Length);
        Route_best(iter,:) = Table(min_index,:);
    else
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter-1), min_Length);
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter,:) = Table(min_index,:);
        else
            Route_best(iter,:) = Route_best(iter-1,:);
        end
    end

    % 更新信息素
    Delta_Tau = zeros(n,n);
    for i = 1:m
        for j = 1:(n-1)
            Delta_Tau(Table(i,j), Table(i,j+1)) = Delta_Tau(Table(i,j), Table(i,j+1)) + Q / Length(i);
        end
        Delta_Tau(Table(i,n), Table(i,1)) = Delta_Tau(Table(i,n), Table(i,1)) + Q / Length(i);
    end
    Tau = (1 - rho) * Tau + Delta_Tau;

    % 迭代次数加一，清空路径记录表
    iter = iter + 1;
    Table = zeros(m, n);
end

%% 结果显示
[Shortest_Length, index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离：', num2str(Shortest_Length)]);
disp(['最短路径：', num2str([Shortest_Route Shortest_Route(1)])]);

%% 绘图
figure(1)
plot([citys(Shortest_Route,1); citys(Shortest_Route(1),1)], ...
    [citys(Shortest_Route,2); citys(Shortest_Route(1),2)], 'o-');
grid on
for i = 1:size(citys, 1)
    text(citys(i, 1), citys(i, 2), ['  ' num2str(i)]);
end
text(citys(Shortest_Route(1),1), citys(Shortest_Route(1),2), '      起点');
text(citys(Shortest_Route(end),1), citys(Shortest_Route(end),2), '      终点');
xlabel('城市位置横坐标');
ylabel('城市位置纵坐标');
title(['蚁群算法路径(最短距离：' num2str(Shortest_Length) ')']);

figure(2)
plot(1:iter_max, Length_best, 'b', 1:iter_max, Length_ave, 'r');
legend('最短距离', '平均距离');
xlabel('迭代次数');
ylabel('距离');
title('各代最短距离与平均距离对比');
