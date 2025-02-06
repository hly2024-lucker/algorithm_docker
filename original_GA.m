clear
clc
close all

X = load('bier127.txt');  % 加载城市坐标数据
X(:,1) = [];            % 移除第一列（假设是城市编号）

% 33958.5975 //0.02

NIND = 45;             % 种群规模
MAXGEN = 8000;          % 最大代数
Pc = 0.9;               % 交叉概率
Pm = 0.02;              % 变异概率
GGAP = 0.9;             % 选择比例
D = Distance(X);        % 生成距离矩阵
N = size(D,1);          % 获取城市数量（假设为34）

%% 初始化种群
Chrom = InitPop(NIND, N);

%% 绘制随机路径图
DrawPath(Chrom(1,:), X)
pause(0.0001);

%% 输出随机路径及其总距离
disp('初始化种群中的一个随机路径:');
OutputPath(Chrom(1,:));
Rlength = PathLength(D, Chrom(1,:));
disp(['总距离:', num2str(Rlength)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%% 优化过程
gen = 0;
figure;
hold on; box on;
xlim([0, MAXGEN]);
title('优化过程');
xlabel('代数');
ylabel('最优路径长度');
ObjV = PathLength(D, Chrom);  % 计算路径长度
preObjV = min(ObjV);
max_N = MAXGEN / 8;  % 连续保持最低进化率的最大代数
total = 1;
tmp = 0;
while gen < MAXGEN
    %% 计算适应度
    ObjV = PathLength(D, Chrom);  % 计算路径长度
    line([gen-1, gen], [preObjV, min(ObjV)]);
    pause(0.0001);
    preObjV = min(ObjV);
    % 新增保留上代最高适应度
    preFitnV = tmp; 
    FitnV = Fitness(ObjV);
    tmp = max(FitnV);
    if tmp <= preFitnV
        total = total + 1;
        fprintf('+1\n');
        if total == max_N
            fprintf("进化停滞，执行下一算法\n");
            break;
        end
    else
        total = 1;
    end

    %% 选择操作
    SelCh = Select(Chrom, FitnV, GGAP);
    % 交叉操作
    SelCh = original_Recombin(SelCh, Pc);
    %% 变异操作
    SelCh = original_Mutate(SelCh, Pm);
    % %% 逆转操作 (可选)
    % SelCh = Reverse(SelCh, D);
    %% 重插入子代到新种群
    Chrom = Reins(Chrom, SelCh, ObjV);
    %% 更新代数
    gen = gen + 1;
end

%% 绘制最优路径图
ObjV = PathLength(D, Chrom);  % 计算路径长度
[minObjV, minInd] = min(ObjV);
Shortest_Route = Chrom(minInd(1), :);

%% 绘图
figure(3);
plot([X(Shortest_Route, 1); X(Shortest_Route(1), 1)], ...
    [X(Shortest_Route, 2); X(Shortest_Route(1), 2)], 'o-');
grid on;
for i = 1:size(X, 1)
    text(X(i, 1), X(i, 2), ['  ' num2str(i)]);
end
text(X(Shortest_Route(1), 1), X(Shortest_Route(1), 2), '      起点');
text(X(Shortest_Route(end), 1), X(Shortest_Route(end), 2), '      终点');
xlabel('城市位置横坐标');
ylabel('城市位置纵坐标');
title(['基于传统遗传算法的路径优化(最短距离：' num2str(minObjV) ')']);

%% 输出最优路径及总距离
disp('最优路径：');
p = OutputPath(Chrom(minInd(1), :));
disp(['总距离:', num2str(ObjV(minInd(1)))]); 
disp('-----------------------------------------------');
