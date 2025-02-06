function [Leader_score,Leader_pos,Convergence_curve]=mu_lambda_BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost)

disp(pop(3).Cost);
%% Problem Definition  �Ż����ⶨ��

CostFunction=@(x) ComputeCost(x);        % Cost Function

%nVar=5;             % Number of Decision Variables  ���߱�������Ŀ�����¡�ʪ�ȵ���Ŀ

VarSize=[1 nVar];   % Decision Variables Matrix Size  ���߾���
VarMin = varMin(1);
VarMax = varMax(1);

%VarMin=-10;         % Decision Variables Lower Bound  ���߱������½�
%VarMax= 10;         % Decision Variables Upper Bound  ���߱������Ͻ�

%% Initialization

% Empty Habitat
habitat.Position=[];
habitat.Cost=[];

% Create Habitats Array
%%repmat(habitat, nPop, 1) ������һ������ nPop ����Ϣ�أ�habitat �ṹ��ĸ��ƣ������飬
%%ÿ����Ϣ�ض����� Position �� Cost �����ֶ�
pop=repmat(habitat,nPop,1);

% Initialize Habitats
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize); %��ʼ��������Ϣ�صĻ���
    pop(i).Cost=CostFunction(pop(i).Position);
end


%% BBO Parameters    BBO�㷨����

%MaxIt=600;          % Maximum Number of Iterations ����������


KeepRate=0.2;                   % Keep Rate    ��������Ϣ�ر���
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats    ÿ����������Ϣ������/��������

nNew=nPop-nKeep;                % Number of New Habitats    ����Ϣ����Ŀ
same_iter = 0;

% Sort Population
[~, SortOrder]=sort([pop.Cost]);
pop=pop(SortOrder);

fitnessValues = [pop.Cost]; 
% % Migration Rates
% mu = B * sin(omega * fitnessValues +varphi)+k2;          % Emmigration Rates    Ǩ���ʣ��������Ǻ���ģ�ͼ���
% lambda = A * cos(omega * fitnessValues +varphi)+k1;  % Immigration Rates    Ǩ���ʣ��������Ǻ���ģ�ͼ���

%% ����ģ�Ͳ�������
a = 1; % Ǩ��������ģ�Ͳ���a���ɸ���ʵ���������
b = 1; % Ǩ��������ģ�Ͳ���b���ɸ���ʵ���������
c = 1; % Ǩ��������ģ�Ͳ���c���ɸ���ʵ���������
A = 0;
B = 0.3;
omega = 1;
varphi = 0;
k1 = 0;
k2 = 0.3;

% % ����ÿ����Ϣ�ص���Ӧ���������ֵ���������Costֵ���棬�ɸ���ʵ���Ż������׼ȷ���壩
% fitnessValues = [pop.Cost]; 
% 
% % Migration Rates
% mu = c * fitnessValues;          % Emmigration Rates    Ǩ���ʣ���������ģ�ͼ���
% lambda = a - b * fitnessValues;  % Immigration Rates    Ǩ���ʣ���������ģ�ͼ���


alpha=0.9;  %Ǩ������,������Ϣ��֮���Ǩ�Ƴ̶�

pMutation=0.02;  %ͻ�����
sigma=0.02*(VarMax-VarMin); %ͻ�����
probability = 0.5;  %Ⱥ�����ƶ���ֵ

%% Migration Rates ԭʼǨ��������
mu=linspace(1,0,nPop);          % Emmigration Rates    Ǩ����
%%����һ������Ϊ nPop �������������� 1 �� 0 �ȼ������ֵ��
%%mu ����Ǩ���ʣ�Emmigration rates��������ʾÿ����Ϣ�ص���Ǩ���ʣ�
%%�� BBO �㷨�У���Ǩ�������𽥼��ٵģ��������������ֵ�Ǵ� 1 �𽥼�С�� 0��

lambda=1-mu;                    % Immigration Rates    Ǩ����


%%

% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop
AvgCost = zeros(MaxIt, 1);  % ����һ�����������洢ÿһ����ƽ������

for it=1:MaxIt
    
    if(same_iter>MaxIt/50)
        disp('_ _ _ _ _ __ __ __ _ _ _ - - -- -  - -- - - ')
        fitnessValues = [pop.Cost]; 
        %%�������ֲ����Ž⣬ת������Ǩ��ģ��
        maxfit = max(fitnessValues);
        minfit = min(fitnessValues);
        fitness_chance = (maxfit - fitnessValues)/(maxfit - minfit);
        
        mu = 0.7*(a - b * fitness_chance)+0.3*(B * sin(omega * (1-fitness_chance) +varphi)+k2);          % Emmigration Rates    Ǩ���ʣ���������ģ�ͼ���
        lambda = 0.7*(c * fitness_chance)+0.3*(A * cos(omega * (1-fitness_chance) +varphi)+k1);  % Immigration Rates    Ǩ���ʣ���������ģ�ͼ���
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
    
    
    
    newpop=pop;
    
    %%�������ƶȲ�����Ƿ���Ҫ�����������ʣ�
%     similarityCount = SimilarityCalculation(pop,nPop);
%     if similarityCount > nPop * 0.4  % �����Ⱥ�г��� 40% �ĸ��������
%         disp('High similarity detected, applying diversity adjustment...');
%         % ��������������������ԣ�������ͻ����ʻ�����Ǩ��
%         pMutation = 0.055;  % ����ͻ�����
%     else
%         pMutation = 0.02;  % �ָ�������ͻ�����
%     end
    for i=1:nPop
        for k=1:nVar
            % Migration  ��⵱ǰ��Ϣ���Ƿ����Ǩ��
            if rand<=lambda(i)
                % Emmigration Probabilities
                EP=mu;
                EP(i)=0;    %�ų���ǰ��Ϣ��
                EP=EP/sum(EP);  %����Ǩ�����������Ϣ���еĸ���
                % Select Source Habitat
                j=RouletteWheelSelection(EP);   %ѡ��Դ��Ϣ�ؽ���Ǩ��
                
                % Migration ���µ�ǰλ��
                %pop(j).Position(k) - pop(i).Position(k) ��ʾԴ��Ϣ�� j �뵱ǰ��Ϣ�� i ֮��Ĳ��졣
                %%�����ϣ���һ�д��뽫��ǰ��Ϣ����Դ��Ϣ��Ǩ�ƣ�����������ȫ���棬
                %%�����ܵ� alpha ��Ӱ�죬ʹ�õ�ǰ��Ϣ�ص�λ�ø��¡�
                newpop(i).Position(k)=pop(i).Position(k) ...
                    +alpha*(pop(j).Position(k)-pop(i).Position(k));

                
            end
    
            % Mutation  ��ԭλ�������һ������Ŷ�
            if rand<=pMutation
                t =max(VarMin,min(VarMax,newpop(i).Position(k)+sigma*randn));
                newpop(i).Position(k) = t;
            end
        end
        
        % Apply Lower and Upper Bound Limits �������½磨�������ޣ�С�����ޣ�
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);
        
        % Evaluation  %������λ�õĳɱ�
        newpop(i).Cost=CostFunction(newpop(i).Position);
    end
    
    % Sort New Population
    [~, SortOrder]=sort([newpop.Cost]);
    newpop=newpop(SortOrder);
    
    % Select Next Iteration Population  ��Ӣ����
    pop=[pop(1:nKeep);
         newpop(1:nNew)];
     
    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Update Best Solution Ever Found
    BestSol=pop(1);
    
    %%�������ƴ����ĸ��ͨ���ж����������ֲ�ͣ��������������
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    if it > 1 && abs(BestCost(it-1) - BestSol.Cost) < 1.0*BestCost(it-1)/MaxIt
        same_iter = same_iter+1;
%         if mod(same_iter,10) == 0
        if same_iter>10
            pMutation = 0.055;  % ����ͻ�����
            sigma = 0.03 * (VarMax - VarMin);  % ���ӱ��췶Χ
        end
    else
        same_iter = 0;
    end
    
    
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Calculate Average Cost for the current generation
    avgCost = mean([pop.Cost]);
    AvgCost(it) = avgCost;  % ����ÿһ�����и����ƽ������

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ', Average Cost = ' num2str(avgCost)]);
    
end

%% Results
% %% Results
% figure;
% 
% % ������Ѵ��۵�ͼ
% semilogy(BestCost, 'LineWidth', 2);
% hold on;  % ���ֵ�ǰͼ�Σ��Ա���ͬһͼ�ϻ�����������
% 
% % ����ƽ�����۵�ͼ
% semilogy(AvgCost, 'LineWidth', 2, 'LineStyle', '-');  % ʹ�����ߣ�'--'����������������
% 
% % ���ñ�ǩ
% xlabel('Iteration');
% ylabel('Cost');
% grid on;  % ������
% 
% % ���ͼ��
% legend('Best Cost', 'Average Cost');

Leader_pos =  BestSol.Position;
Leader_score =  BestSol.Cost;
Convergence_curve = BestCost;
end

