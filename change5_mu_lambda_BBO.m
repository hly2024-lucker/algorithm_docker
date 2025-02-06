function [Leader_score,Leader_pos,Convergence_curve]=change5_mu_lambda_BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost1,ComputeCost2,ComputeCost3)

%% Problem Definition  优化问题定义

CostFunction1=@(x) ComputeCost1(x);        % Cost Function
CostFunction2=@(x) ComputeCost2(x);        % Cost Function
CostFunction3=@(x) ComputeCost3(x);        % Cost Function

Max_Cost1 = 0;  % 存储函数1最大代价
Max_Cost2 = 0;  % 存储函数2最大代价
Max_Cost3 = 0;  % 存储函数3最大代价

Min_Cost1 = pop(1).Cost1;
Min_Cost2 = pop(1).Cost2;
Min_Cost3 = pop(1).Cost3;


VarSize=[1 nVar];   % Decision Variables Matrix Size  决策矩阵
VarMin = varMin(1);
VarMax = varMax(1);

%VarMin=-10;         % Decision Variables Lower Bound  决策变量的下界
%VarMax= 10;         % Decision Variables Upper Bound  决策变量的上界





%% BBO Parameters    BBO算法参数

%MaxIt=600;          % Maximum Number of Iterations 最大迭代次数


KeepRate=0.2;                   % Keep Rate    保留的栖息地比例
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats    每代保留的栖息地数量/四舍五入
nNew=nPop-nKeep;                % Number of New Habitats    新栖息地数目
same_iter = 0;


% Initialize Ideal Points
PosIdealPoint = inf(1, nVar);  % 正理想点
NegIdealPoint = -inf(1, nVar); % 负理想点

fitnessValues = [pop.Cost]; 
% % Migration Rates
% mu = B * sin(omega * fitnessValues +varphi)+k2;          % Emmigration Rates    迁出率，根据三角函数模型计算
% lambda = A * cos(omega * fitnessValues +varphi)+k1;  % Immigration Rates    迁入率，根据三角函数模型计算

%% 线性模型参数设置
a = 1; % 迁入率线性模型参数a，可根据实际情况调整
b = 1; % 迁入率线性模型参数b，可根据实际情况调整
c = 1; % 迁出率线性模型参数c，可根据实际情况调整
A = 0;
B = 0.3;
omega = 1;
varphi = 0;
k1 = 0;
k2 = 0.3;

% % 计算每个栖息地的适应度相关属性值（这里简单用Cost值代替，可根据实际优化问题更准确定义）
% fitnessValues = [pop.Cost]; 
% 
% % Migration Rates
% mu = c * fitnessValues;          % Emmigration Rates    迁出率，根据线性模型计算
% lambda = a - b * fitnessValues;  % Immigration Rates    迁入率，根据线性模型计算


alpha=0.9;  %迁移因子,控制栖息地之间的迁移程度

% pMutation=0.02;  %突变概率
% sigma=0.02*(VarMax-VarMin); %突变幅度
probability = 0.85;  %群体相似度阈值
pMutationInitial = 0.02;  %突变概率
sigma=0.02*(VarMax-VarMin); %突变幅度

%三种函数在总代价函数中权重
m1 = 0.3;
m2 = 0.4;
m3 = 0.3;




%% Migration Rates 原始迁移率设置
mu=linspace(1,0,nPop);          % Emmigration Rates    迁出率
%%生成一个长度为 nPop 的向量，包含从 1 到 0 等间隔的数值。
%%mu 是外迁概率（Emmigration rates），它表示每个栖息地的外迁概率，
%%在 BBO 算法中，外迁概率是逐渐减少的，所以这个向量的值是从 1 逐渐减小到 0。

lambda=1-mu;                    % Immigration Rates    迁入率


%%

% % Empty Habitat
% habitat.Position=[];
% habitat.Cost1=[];
% habitat.Cost2=[];
% habitat.Cost3=[];
% habitat.Cost=[];


BestCost=zeros(MaxIt,1);
BestCost1=zeros(MaxIt,1);
BestCost2=zeros(MaxIt,1);
BestCost3=zeros(MaxIt,1);
% 
% % Create Habitats Array
% %%repmat(habitat, nPop, 1) 会生成一个包含 nPop 个栖息地（habitat 结构体的复制）的数组，
% %%每个栖息地都包含 Position 和 Cost 两个字段
% pop=repmat(habitat,nPop,1);
% Initialize Habitats
for i=1:nPop
%     pop(i).Position=unifrnd(VarMin,VarMax,VarSize); %初始化各个栖息地的环境
%     
%     pop(i).Cost1=CostFunction1(pop(i).Position);
%     pop(i).Cost2=CostFunction2(pop(i).Position);
%     pop(i).Cost3=CostFunction3(pop(i).Position);
    if pop(i).Cost1>Max_Cost1
        Max_Cost1 = pop(i).Cost1;
    end
    if pop(i).Cost2>Max_Cost2
        Max_Cost2 = pop(i).Cost2;
    end
    if pop(i).Cost3>Max_Cost3
        Max_Cost3 = pop(i).Cost3;
    end
    if pop(i).Cost1<Min_Cost1
        Min_Cost1 = pop(i).Cost1;
    end
    if pop(i).Cost2< Min_Cost2
        Min_Cost2 = pop(i).Cost2;
    end
    if pop(i).Cost3<Min_Cost3
        Min_Cost3 = pop(i).Cost3;
    end
    
%     disp(['cost1 is :' num2str(Max_Cost1)  'cost2 is :' num2str(Max_Cost2)  'cost3 is :' num2str(Max_Cost3)]);
end

for i=1:nPop
    %初始化计算种群内个体的总代价
    pop(i).Cost = m1*(Max_Cost1-Min_Cost1)/(Max_Cost1-pop(i).Cost1) + m2*(Max_Cost2-Min_Cost2)/(Max_Cost2-pop(i).Cost2) + m3*(Max_Cost3-Min_Cost3)/(Max_Cost3-pop(i).Cost3);
end

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop
AvgCost = zeros(MaxIt, 1);  % 创建一个空数组来存储每一代的平均代价

for it=1:MaxIt
    for i=1:nPop
        for j=1:nVar
            
            %% 动态更新正负理想点 （理想点的每个变量的值都是最小/最大的）
            if abs(pop(i).Position(j)) < abs(PosIdealPoint(j))
                PosIdealPoint(j) = pop(i).Position(j);
            end
            if abs(pop(i).Position(j)) > abs(NegIdealPoint(j))
                NegIdealPoint(j) = pop(i).Position(j);
            end
        end
    end
    if(same_iter>MaxIt/50)
        disp('_ _ _ _ _陷入局部最优解，转入新型迁移模型___ _   __  ____  ___ ')
%        fitnessValues = [pop.Cost]; 
        
        for i = 1:nPop
            distanceToPosIdeal = sqrt(sum((abs(pop(i).Position) - abs(PosIdealPoint)).^2));
            distanceToNegIdeal = sqrt(sum((abs(pop(i).Position) - abs(NegIdealPoint)).^2));
            fitnessValues(i) = distanceToPosIdeal / (distanceToNegIdeal + eps); % 更新适应度
        end
    
        
        %%如果陷入局部最优解，转入新型迁移模型
        maxfit = max(fitnessValues);
        minfit = min(fitnessValues);
        fitness_chance = (maxfit - fitnessValues)/(maxfit - minfit);
        
        mu = 0.7*(a - b * fitness_chance)+0.3*(B * sin(omega * (1-fitness_chance) +varphi)+k2);          % Emmigration Rates    迁出率，根据线性模型计算
        lambda = 0.7*(c * fitness_chance)+0.3*(A * cos(omega * (1-fitness_chance) +varphi)+k1);  % Immigration Rates    迁入率，根据线性模型计算
%        lambda = 1 - mu; 
    
        min_mu = min(mu);
        max_mu = max(mu);
        mu_mapped = (mu - min_mu) / (max_mu - min_mu);
        
        min_lambda = min(lambda);
        max_lambda = max(lambda);
        lambda_mapped = (lambda - min_lambda) / (max_lambda - min_lambda);
        
        mu = mu_mapped;
        lambda = lambda_mapped;
        
    end
    
    
    %%设置每个栖息地的变异概率
    costThreshold = median([pop.Cost]); % 成本阈值，这里取中位数
    for i = 1:length(pop)
        if pop(i).Cost > costThreshold
            pop(i).pMutation = pMutationInitial * 1.3; % 适应度差的栖息地变异率加da
        else
            pop(i).pMutation = pMutationInitial * 0.9;
        end
    end
    
    newpop=pop;
    
    %%去重操作
%     newpop = Remove_Calculation(pop,nPop,VarMin,VarMax,CostFunction);
Remove_Threshold = 0.1;  %  去重阈值

for i = 1:nPop
    for j = i+1:nPop
        % 计算个体 i 和 j 之间的相似度
        dist = DistanceCalculation(newpop(i).Position, newpop(j).Position);
        % 如果距离小于某个阈值，则认为它们相似
        if dist < Remove_Threshold
            p1 = local_search(CostFunction1,newpop(1).Position,VarMin,VarMax);
            p2 = local_search(CostFunction2,newpop(1).Position,VarMin,VarMax);
            p3 = local_search(CostFunction3,newpop(1).Position,VarMin,VarMax);
            
            
            numCost1=CostFunction1(p1);
            numCost2=CostFunction2(p1);
            numCost3=CostFunction3(p1);
            
            endCost1 = m1*(Max_Cost1-Min_Cost1)/(Max_Cost1-numCost1) + m2*(Max_Cost2-Min_Cost2)/(Max_Cost2-numCost2) + m3*(Max_Cost3-Min_Cost3)/(Max_Cost3-numCost3);
            
            numCost1=CostFunction1(p2);
            numCost2=CostFunction2(p2);
            numCost3=CostFunction3(p2);
            
            endCost2 = m1*(Max_Cost1-Min_Cost1)/(Max_Cost1-numCost1) + m2*(Max_Cost2-Min_Cost2)/(Max_Cost2-numCost2) + m3*(Max_Cost3-Min_Cost3)/(Max_Cost3-numCost3);

            numCost1=CostFunction1(p3);
            numCost2=CostFunction2(p3);
            numCost3=CostFunction3(p3);
            
            endCost3 = m1*(Max_Cost1-Min_Cost1)/(Max_Cost1-numCost1) + m2*(Max_Cost2-Min_Cost2)/(Max_Cost2-numCost2) + m3*(Max_Cost3-Min_Cost3)/(Max_Cost3-numCost3);
            tcost = endCost1;
            tindex = 1;
            if tcost < endCost2
                tcost = endCost2;
                tindex = 2;
            end
            if tcost < endCost3
                tcost = endCost3;
                tindex = 3;
            end
            
            if tcost<newpop(i).Cost
                if tindex == 1
                    newpop(j).Position = p1;
                elseif tindex==2
                    newpop(j).Position = p2;
                else
                    newpop(j).Position = p3;
                end
            newpop(j).Cost = tcost;      
            end
            
    if newpop(j).Cost1>Max_Cost1
        Max_Cost1 = newpop(i).Cost1;
    end
    if newpop(j).Cost2>Max_Cost2
        Max_Cost2 = newpop(i).Cost2;
    end
    if newpop(j).Cost3>Max_Cost3
        Max_Cost3 = newpop(i).Cost3;
    end
    if pop(i).Cost1<Min_Cost1
        Min_Cost1 = pop(i).Cost1;
    end
    if pop(i).Cost2< Min_Cost2
        Min_Cost2 = pop(i).Cost2;
    end
    if pop(i).Cost3<Min_Cost3
        Min_Cost3 = pop(i).Cost3;
    end
    
% %             pcost = ComputeCost(p);
%             %如果局部搜索后的解更优，则取代最优解
%             pop(j).Position = p;
%             pop(j).Cost = pcost;
        end
    end
end    


    
    %%计算相似度并检查是否需要调整（变异率）
    similarityCount = SimilarityCalculation(pop,nPop);
    if similarityCount > nPop * probability  % 如果种群中超过 40% 的个体对相似
        disp('High similarity detected, applying diversity adjustment...');
        % 可以在这里引入调整策略，如增加突变概率或引导迁移
        pMutationInitial = 0.055;  % 增加突变概率
    else
        pMutationInitial = 0.02+0.35*(it/MaxIt);  % 恢复正常的突变概率
    end
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
            
%             % Mutation  在原位置上添加一个随机扰动
%             if rand<=pMutation
%                 t =max(VarMin,min(VarMax,newpop(i).Position(k)+sigma*randn));
%                 newpop(i).Position(k) = t;
%             end
            % Mutation  在原位置上添加一个随机扰动
            if rand<=pop(i).pMutation
                t =max(VarMin,min(VarMax,newpop(i).Position(k)+sigma*randn));
                newpop(i).Position(k) = t;
            end

        end
        
        % Apply Lower and Upper Bound Limits 限制上下界（大于下限，小于上限）
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
%         % Evaluation  %计算新位置的成本
%         newpop(i).Cost=CostFunction(newpop(i).Position);
        newpop(i).Cost1=CostFunction1(newpop(i).Position);
        newpop(i).Cost2=CostFunction2(newpop(i).Position);
        newpop(i).Cost3=CostFunction3(newpop(i).Position);
        if newpop(i).Cost1>Max_Cost1
            Max_Cost1 = newpop(i).Cost1;
        end
        if newpop(i).Cost2>Max_Cost2
            Max_Cost2 = newpop(i).Cost2;
        end
        if newpop(i).Cost3>Max_Cost3
            Max_Cost3 = newpop(i).Cost3;
        end
        if pop(i).Cost1<Min_Cost1
            Min_Cost1 = pop(i).Cost1;
        end
        if pop(i).Cost2< Min_Cost2
            Min_Cost2 = pop(i).Cost2;
        end
        if pop(i).Cost3<Min_Cost3
            Min_Cost3 = pop(i).Cost3;
        end
        
        newpop(i).Cost = m1*(Max_Cost1-Min_Cost1)/(Max_Cost1-newpop(i).Cost1) + m2*(Max_Cost2-Min_Cost2)/(Max_Cost2-newpop(i).Cost2) + m3*(Max_Cost3-Min_Cost3)/(Max_Cost3-newpop(i).Cost3);
    end
    
    disp('---------------');
    %%利用理想点计算适应度
    for i = 1:nPop
        distanceToPosIdeal = sqrt(sum((abs(newpop(i).Position) - abs(PosIdealPoint)).^2));
        distanceToNegIdeal = sqrt(sum((abs(newpop(i).Position) - abs(NegIdealPoint)).^2));
        fitnessValues(i) = distanceToPosIdeal / (distanceToNegIdeal + eps); % 更新适应度
    end

    [~, SortOrder]=sort(fitnessValues);     %这里以与正负理想点的距离远近程度来排序
    newpop=newpop(SortOrder);
    fitnessValues = fitnessValues(SortOrder);
    disp(['ow :' num2str(newpop(1).Cost) '-- ' num2str(newpop(1).Cost1) '+' num2str(newpop(i).Cost2) '+' num2str(newpop(i).Cost3)]);

    % Sort Population
    [~, SortOrder]=sort([pop.Cost],'descend');
    pop=pop(SortOrder);
    disp(['pow :' num2str(pop(1).Cost) '-- ' num2str(pop(1).Cost1) '+' num2str(pop(i).Cost2) '+' num2str(pop(i).Cost3)]);
    pop=[pop(1:nKeep);
         newpop(1:nNew)];
     
    % Update Best Solution Ever Found
    BestSol=pop(1);
    
    %%引入相似代数的概念，通过判断连续几代局部停滞来调整变异率
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    if it > 1 && abs(BestCost(it-1) - BestSol.Cost) < 1.0*BestCost(it-1)/MaxIt
        same_iter = same_iter+1;
        if same_iter>10
            pMutationInitial = 0.055;  % 增加突变概率
            sigma = 0.03 * (VarMax - VarMin);  % 增加变异范围
        end
    else
        same_iter = 0;
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Calculate Average Cost for the current generation
    avgCost = mean([pop.Cost]);
    AvgCost(it) = avgCost;  % 计算每一代所有个体的平均代价

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': c5Best Cost = ' num2str(BestCost(it)) ', Average Cost = ' num2str(avgCost)]);
    
end

Leader_pos =  BestSol.Position;
Leader_score =  BestSol.Cost;
Convergence_curve = BestCost;
end

