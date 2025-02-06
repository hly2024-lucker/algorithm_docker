function [Leader_score,Leader_pos,Convergence_curve1,Convergence_curve2,Convergence_curve3]=BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost1,ComputeCost2,ComputeCost3)

disp(pop(3).Cost);



%% Problem Definition  优化问题定义

CostFunction1=@(x) ComputeCost1(x);        % Cost Function
CostFunction2=@(x) ComputeCost2(x);        % Cost Function
CostFunction3=@(x) ComputeCost3(x);        % Cost Function
% 
% Max_Cost1 = zeros(MaxIt, 1);  % 创建一个空数组来存储每一代的函数1的最大代价
% Max_Cost2 = zeros(MaxIt, 1);  % 创建一个空数组来存储每一代的函数2最大代价
% Max_Cost3 = zeros(MaxIt, 1);  % 创建一个空数组来存储每一代的函数3最大代价

Max_Cost1 = 0;  % 存储函数1最大代价
Max_Cost2 = 0;  % 存储函数2最大代价
Max_Cost3 = 0;  % 存储函数3最大代价


VarSize=[1 nVar];   % Decision Variables Matrix Size  决策矩阵
VarMin = varMin(1);
VarMax = varMax(1);

%VarMin=-10;         % Decision Variables Lower Bound  决策变量的下界
%VarMax= 10;         % Decision Variables Upper Bound  决策变量的上界

%% BBO Parameters    BBO算法参数

KeepRate=0.2;                   % Keep Rate    保留的栖息地比例
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats    每代保留的栖息地数量/四舍五入
nNew=nPop-nKeep;                % Number of New Habitats    新栖息地数目

% Migration Rates
mu=linspace(1,0,nPop);          % Emmigration Rates    迁出率
%%生成一个长度为 nPop 的向量，包含从 1 到 0 等间隔的数值。
%%mu 是外迁概率（Emmigration rates），它表示每个栖息地的外迁概率，
%%在 BBO 算法中，外迁概率是逐渐减少的，所以这个向量的值是从 1 逐渐减小到 0。
lambda=1-mu;                    % Immigration Rates    迁入率
alpha=0.9;  %迁移因子,控制栖息地之间的迁移程度
pMutation=0.1;  %突变概率
sigma=0.02*(VarMax-VarMin); %突变幅度
probability = 0.5;  %群体相似度阈值

%三种函数在总代价函数中权重
m1 = 0.3;
m2 = 0.4;
m3 = 0.3;

%% Initialization
% 
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


% Create Habitats Array
%%repmat(habitat, nPop, 1) 会生成一个包含 nPop 个栖息地（habitat 结构体的复制）的数组，
%%每个栖息地都包含 Position 和 Cost 两个字段
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
end
for i=1:nPop
    %初始化计算种群内个体的总代价
    pop(i).Cost = m1*pop(i).Cost1 + m2*pop(i).Cost2 + m3*pop(i).Cost3;
end 
%% BBO Main Loop
AvgCost = zeros(MaxIt, 1);  % 创建一个空数组来存储每一代的平均代价

for it=1:MaxIt
    
    newpop=pop;

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
            if rand<=pMutation
                t =max(VarMin,min(VarMax,newpop(i).Position(k)+sigma*randn));
                newpop(i).Position(k) = t;
            end
        end
        
        % Apply Lower and Upper Bound Limits 限制上下界（大于下限，小于上限）
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
        % Evaluation  %计算新位置的成本
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
    end
    
    % 计算每代每个函数下的最小值
    [~, SortOrder]=sort([newpop.Cost1]);
    BestCost1(it) = newpop(1).Cost1;
    [~, SortOrder]=sort([newpop.Cost2]);
    BestCost2(it) = newpop(1).Cost2;
    [~, SortOrder]=sort([newpop.Cost3]);
    BestCost3(it) = newpop(1).Cost3;
    
    for i=1:nPop
        %%计算当前群体的混合代价
        newpop(i).Cost = m1*newpop(i).Cost1 + m2*newpop(i).Cost2 + m3*newpop(i).Cost3;
    end
    disp([ 'maxcost1 is: ' num2str(Max_Cost1) ', maxcost2 is: ' num2str(Max_Cost2) ',maxcost3 is: ' num2str(Max_Cost3)]);
    [~, SortOrder]=sort([newpop.Cost],'descend');
    newpop=newpop(SortOrder);
    disp(['ow :' num2str(newpop(1).Cost) '-- ' num2str(newpop(1).Cost1) '+' num2str(newpop(i).Cost2) '+' num2str(newpop(i).Cost3)]);
    
    

    % Sort Population
    [~, SortOrder]=sort([pop.Cost],'descend');
    pop=pop(SortOrder);
    disp(['pow :' num2str(pop(1).Cost) '-- ' num2str(pop(1).Cost1) '+' num2str(pop(i).Cost2) '+' num2str(pop(i).Cost3)]);

    % Select Next Iteration Population  精英策略
    pop=[pop(1:nKeep);
         newpop(1:nNew)];
     
    [~, SortOrder]=sort([pop.Cost],'descend');
    pop=pop(SortOrder);

    % Update Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BestCost1(it) = pop(1).Cost1;
    BestCost2(it) = pop(1).Cost2;
    BestCost3(it) = pop(1).Cost3;

    % Calculate Average Cost for the current generation
    avgCost = mean([pop.Cost]);
    AvgCost(it) = avgCost;  % 计算每一代所有个体的平均代价
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': bboBest Cost = ' num2str(BestCost(it)) ', Average Cost = ' num2str(avgCost)]);
    
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
Convergence_curve1 = BestCost1;
Convergence_curve2 = BestCost2;
Convergence_curve3 = BestCost3;

end

