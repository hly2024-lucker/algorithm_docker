function [Leader_score,Leader_pos,Convergence_curve]=BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost)

%%%%%%只修改单个栖息地变异概率（较差的部分变异概率大），迁移率跟原始无变化
%% Problem Definition  优化问题定义

CostFunction=@(x) ComputeCost(x);        % Cost Function

%nVar=5;             % Number of Decision Variables  决策变量的数目：气温、湿度等数目

VarSize=[1 nVar];   % Decision Variables Matrix Size  决策矩阵
VarMin = varMin(1);
VarMax = varMax(1);

%VarMin=-10;         % Decision Variables Lower Bound  决策变量的下界
%VarMax= 10;         % Decision Variables Upper Bound  决策变量的上界
%%按照适应度对初始种群排序，第一个最好
[sortedCost, sortOrder] = sort([pop.Cost]);
pop = pop(sortOrder);


%% BBO Parameters    BBO算法参数

KeepRate=0.2;                   % Keep Rate    保留的栖息地比例
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats    每代保留的栖息地数量/四舍五入

nNew=nPop-nKeep;                % Number of New Habitats    新栖息地数目

% 迁入迁出率设置
mu=linspace(1,0,nPop);          % Emmigration Rates    迁出率
%%生成一个长度为 nPop 的向量，包含从 1 到 0 等间隔的数值。
%%mu 是外迁概率（Emmigration rates），它表示每个栖息地的外迁概率，
%%在 BBO 算法中，外迁概率是逐渐减少的，所以这个向量的值是从 1 逐渐减小到 0。
lambda=1-mu;                    % Immigration Rates    迁入率

%%原始参数设置
alpha=0.9;  %迁移因子,控制栖息地之间的迁移程度
lambdaFactor = 0.9; % 调整因子，可根据实验调整
muFactor = 1.1; % 调整因子，可根据实验调整
pMutationInitial = 0.02;  %突变概率
sigma=0.02*(VarMax-VarMin); %突变幅度
probability = 0.5;  %群体相似度阈值
same_iter = 0;

%%==========================================================================


% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop
AvgCost = zeros(MaxIt, 1);  % 创建一个空数组来存储每一代的平均代价

for it=1:MaxIt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%设置每个栖息地的变异概率
    costThreshold = median([pop.Cost]); % 成本阈值，这里取中位数
    for i = 1:length(pop)
        if pop(i).Cost > costThreshold
            pop(i).pMutation = pMutationInitial * 1.3; % 适应度差的栖息地变异率加da
        else
            pop(i).pMutation = pMutationInitial * 0.7;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newpop=pop;
    
    %%计算相似度并检查是否需要调整（变异率）
    similarityCount = SimilarityCalculation(pop,nPop);
    if similarityCount > nPop * 0.4  % 如果种群中超过 40% 的个体对相似
        disp('High similarity detected, applying diversity adjustment...');
        % 可以在这里引入调整策略，如增加突变概率或引导迁移
        pMutationInitial = 0.055;  % 增加突变概率
        sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围

        mu = mu * 0.8;    % 迁出率
        lambda = lambda * 1.2;
    else
        pMutationInitial = 0.02;  % 恢复正常的突变概率
        sigma = 0.02 * (VarMax - VarMin);

    end
    EP = zeros(1, nPop);
    for i=1:nPop
        for k=1:nVar
            % Migration  检测当前栖息地是否进行迁移
            if rand<=lambda(i)
                % Emmigration Probabilities
                
                EP=mu;
                EP(i)=0;    %排除当前栖息地
                EP=EP/sum(EP);  %计算迁入其余各个栖息地中的概率
                % Select Source Habitat
                j=RouletteWheelSelection(EP);   %选择源栖息地进行迁移
                
                % Migration 更新当前位置
                %pop(j).Position(k) - pop(i).Position(k) 表示源栖息地 j 与当前栖息地 i 之间的差异。
                %%整体上，这一行代码将当前栖息地向源栖息地迁移，但并不是完全跟随，
                %%而是受到 alpha 的影响，使得当前栖息地的位置更新。
                newpop(i).Position(k)=pop(i).Position(k) ...
                    +alpha*(pop(j).Position(k)-pop(i).Position(k));

                
            end
            
            % Mutation  在原位置上添加一个随机扰动
            if rand<=pop(i).pMutation
                t =max(VarMin,min(VarMax,newpop(i).Position(k)+sigma*randn));
                newpop(i).Position(k) = t;
            end
        end
        
        % Apply Lower and Upper Bound Limits 限制上下界（大于下限，小于上限）
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
        % Evaluation  %计算新位置的成本
        newpop(i).Cost=CostFunction(newpop(i).Position);
    end
    
    % Sort New Population
    [~, SortOrder]=sort([newpop.Cost]);
    newpop=newpop(SortOrder);
    
    % Select Next Iteration Population  精英策略
    pop=[pop(1:nKeep);
         newpop(1:nNew)];
     
    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Update Best Solution Ever Found
    BestSol=pop(1);
    %%引入相似代数的概念，通过判断连续几代局部停滞来调整变异率
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    if it > 1 && abs(BestCost(it-1) - BestSol.Cost)==0
        same_iter = same_iter+1;
        disp('_ _ _ _ _ __ __ __ _ _ _ - - -- -  - -- - - ')
        if mod(same_iter,10) == 0
            pMutationInitial = 0.055;  % 增加突变概率
            sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围
        end
    else
        same_iter = 0;
    end
    %%如果连续较多代数相似，执行局部搜索算法，对最优解进行可能的替换  
    %%减小变异概率，使种群快速收敛，从而进行相似度的判定调整
    if same_iter>50
        pMutationInitial = 0.03;  % 增加突变概率
        sigma = 0.01 * (VarMax - VarMin);  % 增加变异范围
        p = local_search(CostFunction,pop(1).Position,VarMin,VarMax);
        pcost = ComputeCost(p);
        %如果局部搜索后的解更优，则取代最优解
        if(pcost<pop(1).Cost)
            pop(1).Position = p;
            pop(1).Cost = pcost;
        end
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

    
    % Calculate Average Cost for the current generation
    avgCost = mean([pop.Cost]);
    AvgCost(it) = avgCost;  % 计算每一代所有个体的平均代价
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ', Average Cost = ' num2str(avgCost)]);
    
end

%% Results
% %% Results
% figure;
% 
% % 绘制最佳代价的图
% semilogy(BestCost, 'LineWidth', 2);
% hold on;  % 保持当前图形，以便在同一图上绘制其他数据
% 
% % 绘制平均代价的图
% semilogy(AvgCost, 'LineWidth', 2, 'LineStyle', '-');  % 使用虚线（'--'）来区分两条曲线
% 
% % 设置标签
% xlabel('Iteration');
% ylabel('Cost');
% grid on;  % 打开网格
% 
% % 添加图例
% legend('Best Cost', 'Average Cost');

Leader_pos =  BestSol.Position;
Leader_score =  BestSol.Cost;
Convergence_curve = BestCost;
end

