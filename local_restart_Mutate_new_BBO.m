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
function [Leader_score,Leader_pos,Convergence_curve]=local_restart_Mutate_new_BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost)
%%[Best_score,Best_pos,cg_curve]

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

pmu = linspace(1, 0, nPop);             % 迁出率
plambda = 1 - pmu;                       % 迁入率

alpha = 0.9;                           % 迁移因子
pMutation = 0.1;                       % 初始变异概率
sigma = 0.02 * (VarMax - VarMin);      % 变异幅度
convergenceThreshold = 1e-400;           % 收敛阈值


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
    %% 自适应迁移率调整(目前已经移除此功能)
    % 根据当前代数调整迁移率
    pmu = linspace(1, 0, nPop)* (1 - (it / MaxIt)*0.8);    % 随代数减小迁出率
    plambda = (1 - pmu)* (1 - (it / MaxIt)*0.8);
    mu = pmu;
    lambda = plambda;
    % 计算种群相似度并调整变异概率
    similarityCount = SimilarityCalculation(pop, nPop);
    if similarityCount > nPop * 0.4
        disp('检测到高相似度，增强种群多样性...');
        pMutation = 0.05 + 0.01 * (1 - it / (MaxIt+1));  %临时增大变异率跳出相似度+随代数增加早期变异强度
        sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围,随代数减小变异幅度
        mu = pmu*0.9;
        lambda = plambda * 1.1;
        if similarityCount > nPop * 0.6
            mu = pmu*0.8;
            lambda = plambda * 1.2;
        end
    else
        pMutation = 0.02 + 0.03 * (1 - it / (MaxIt+1));  % 随代数增加早期变异强度
        sigma = (0.02 + 0.01 * (1 - it / (MaxIt+1)))* (VarMax - VarMin);  % 随代数减小变异幅度
    end

 
    %% 随机重启：当算法停滞时，重新随机选择部分个体  （其中随机重启函数也可以用栖息地初始化的语句pop(i).Position = unifrnd(VarMin, VarMax, VarSize);）
        if(mod(same_iter,MaxIt/3) == 0&&same_iter>0)
            random_restart_idx = randperm(nPop, round(nPop * 0.3));  % 随机选择部分个体  15个
            % 批量赋值Position
            for i = 1:length(random_restart_idx)
                t = random_restart_idx(i);
                pop(t).Position = random_restart(1,nVar,VarMin,VarMax);  % 逐个赋值Position
                pop(t).Cost = ComputeCost(pop(t).Position);
            end
            same_iter = 0;
        end 



    %% 迁移和变异
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

    %% 若上下两代收敛幅度较大，则局部搜索机制进行改进,对最优解局部探索

    if(it>MaxIt*0.9||mod(same_iter,20)==0)
        p = local_search(CostFunction,pop(1).Position,VarMin,VarMax);
        pcost = ComputeCost(p);
        %如果局部搜索后的解更优，则取代最优解
        if(pcost<pop(1).Cost)
            pop(1).Position = p;
            pop(1).Cost = pcost;
        end
    end
    BestSol = pop(1);
    % 存储最佳和平均代价
    BestCost(it) = BestSol.Cost;
    AvgCost(it) = mean([pop.Cost]);
    
    %%%多代不变，进行调整
    tmp = max(it-1,1);
%     if it > 1 && abs(BestCost(tmp) - BestSol.Cost) < abs(BestCost(tmp))*(1/MaxIt)
    if(abs(BestCost(tmp) - BestSol.Cost)==0)
        disp('---------------------------------');
        same_iter = same_iter+1;
        if(same_iter>50)
            pMutation = 0.05 + 0.01 * (1 - it / (MaxIt+1));
            sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围,随代数减小变异幅度
            if(same_iter>100)
                pMutation = 0.05 + 0.03 * (1 - it / (MaxIt+1));
                if(same_iter>150)
                    pMutation = 0.05 + 0.05 * (1 - it / (MaxIt+1));
                end
            end
        end
    else 
        same_iter = 0;
    end

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
   %%[Best_score,Best_pos,cg_curve]  
%[Leader_pos,Leader_score,Convergence_curve] = [BestSol.Cost,BestSol.Position,BestCost];
end
