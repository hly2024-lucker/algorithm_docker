%
% 版权所有 (c) 2015, Yarpiz (www.yarpiz.com)
% 保留所有权利。请阅读“license.txt”以了解许可证条款。
%
% 项目代码: YPEA113
% 项目标题: 改进版 BBO 算法（采用增强的变异策略和收敛阈值）
% 发布者: Yarpiz (www.yarpiz.com)
%
% 开发者: S. Mostapha Kalami Heris (Yarpiz 团队成员)
% 修改者: [您的名字]
%
function [Leader_score,Leader_pos,Convergence_curve]=Mutatenew_BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost)


%% 问题定义
CostFunction = @(x) ComputeCost(x);    % 成本函数

%nVar = 5;                              % 决策变量的数量
VarSize = [1 nVar];                    % 决策变量矩阵的大小
VarMin = varMin(1);                         % 决策变量的下界
VarMax = varMax(1);                           % 决策变量的上界

%% BBO 参数
% MaxIt = 600;                          % 最大迭代次数
% nPop = 50;                             % 种群规模（栖息地数量）

KeepRate = 0.2;                        % 精英保留率
nKeep = round(KeepRate * nPop);        % 每代保留的栖息地数量
nNew = nPop - nKeep;                   % 每代新生成的栖息地数量

mu = linspace(1, 0, nPop);             % 迁出率
lambda = 1 - mu;                       % 迁入率

alpha = 0.9;                           % 迁移因子
pMutation = 0.1;                       % 初始变异概率
sigma = 0.02 * (VarMax - VarMin);      % 变异幅度
convergenceThreshold = 1e-400;           % 收敛阈值

% %% 初始化
% 
% % 空栖息地结构
% habitat.Position = [];
% habitat.Cost = [];
% 
% % 创建初始种群数组
% pop = repmat(habitat, nPop, 1);
% 
% 
% % 初始化栖息地
% for i = 1:nPop
%     pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
%     pop(i).Cost = CostFunction(pop(i).Position);
% end

% 按成本排序种群
[~, SortOrder] = sort([pop.Cost]);
pop = pop(SortOrder);

% 最佳解决方案
BestSol = pop(1);
same_iter = 0;
% 存储最佳代价和平均代价的数组
BestCost = zeros(MaxIt, 1);
AvgCost = zeros(MaxIt, 1);

%% BBO 主循环

for it = 1:MaxIt
    newpop = pop;

    % 计算种群相似度并调整变异概率
    similarityCount = SimilarityCalculation(pop, nPop);
    if similarityCount > nPop * 0.4
        disp('检测到高相似度，增强种群多样性...');
        pMutation = 0.05 + 0.01 * (1 - it / (MaxIt+1));  % 临时增大变异率跳出相似度+随代数增加早期变异强度
        sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围,随代数减小变异幅度
    else
        pMutation = 0.02 + 0.02 * (1 - it / (MaxIt+1));  % 随代数增加早期变异强度
        sigma = 0.02 * (VarMax - VarMin);  % 随代数减小变异幅度
    end

% 基于成本收敛调整变异概率--------收敛性阈值（待改进成连续几代）
%      if it > 1 && abs(BestCost(it-1) - BestSol.Cost) < abs(BestCost(it-1))*(0.00001/MaxIt)
%         disp('-------------------------------');
%         pMutation = min(pMutation*1.01,0.4);       % 若接近收敛则增加变异概率
%     end

    if it > 1 && abs(BestCost(it-1) - BestSol.Cost) < abs(BestCost(it-1))*(0.0001/MaxIt)
        same_iter = same_iter+1;
        disp('______________')
        if mod(same_iter,5) == 0
            pMutation = 0.06;
            sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围,随代数减小变异幅度
        end

    end



%% 自适应迁移率调整
    % 根据当前代数调整迁移率
    mu = linspace(1, 0, nPop) * (1 - (it / MaxIt+1));    % 随代数减小迁移率
    lambda = 1 - mu; 
    
    % 迁移和变异
    for i = 1:nPop
        for k = 1:nVar
            % 迁移操作
            if rand <= lambda(i)
                EP = mu;
                EP(i) = 0;                % 排除当前栖息地
                EP = EP / sum(EP);        % 归一化概率
                % 选择源栖息地并进行迁移
                j = RouletteWheelSelection(EP);
                newpop(i).Position(k) = pop(i).Position(k) + alpha * ( pop(j).Position(k) - pop(i).Position(k) );
          
            end

            % 自适应变异（带有高斯噪声）
            if rand <= pMutation
                newpop(i).Position(k) = newpop(i).Position(k) + sigma * randn;
            end
        end

        % 边界限制
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);

        % 计算新位置的成本
        newpop(i).Cost = CostFunction(newpop(i).Position);
    end

    % 对新种群进行排序
    [~, SortOrder] = sort([newpop.Cost]);
    newpop = newpop(SortOrder);

    % 选择下一代种群
    pop = [pop(1:nKeep); newpop(1:nNew)];

    % 排序种群并更新最佳解决方案
    [~, SortOrder] = sort([pop.Cost]);
    pop = pop(SortOrder);
    BestSol = pop(1);

    % 存储最佳和平均代价
    BestCost(it) = BestSol.Cost;
    AvgCost(it) = mean([pop.Cost]);

    % 显示迭代信息
    disp(['迭代 ' num2str(it) ': 最佳代价 = ' num2str(BestCost(it)) ', 平均代价 = ' num2str(AvgCost(it))]);

    % 检查收敛
    if it > 1 && abs(BestCost(it) - BestCost(it - 1)) < convergenceThreshold
        disp('达到收敛阈值，提前停止...');
        BestCost = BestCost(1:it);
        AvgCost = AvgCost(1:it);
        break;
    end
end

% %% 结果
% figure;
% semilogy(BestCost, 'LineWidth', 2);
% hold on;
% semilogy(AvgCost, 'LineWidth', 2, 'LineStyle', '-');
% xlabel('迭代次数');
% ylabel('代价');
% grid on;
% legend('最佳代价', '平均代价');

Leader_pos =  BestSol.Position;
Leader_score =  BestSol.Cost;
Convergence_curve = BestCost;
end