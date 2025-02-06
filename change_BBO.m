function [Leader_score,Leader_pos,Convergence_curve]=BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost)

%%%%%%ֻ�޸ĵ�����Ϣ�ر�����ʣ��ϲ�Ĳ��ֱ�����ʴ󣩣�Ǩ���ʸ�ԭʼ�ޱ仯
%% Problem Definition  �Ż����ⶨ��

CostFunction=@(x) ComputeCost(x);        % Cost Function

%nVar=5;             % Number of Decision Variables  ���߱�������Ŀ�����¡�ʪ�ȵ���Ŀ

VarSize=[1 nVar];   % Decision Variables Matrix Size  ���߾���
VarMin = varMin(1);
VarMax = varMax(1);

%VarMin=-10;         % Decision Variables Lower Bound  ���߱������½�
%VarMax= 10;         % Decision Variables Upper Bound  ���߱������Ͻ�
%%������Ӧ�ȶԳ�ʼ��Ⱥ���򣬵�һ�����
[sortedCost, sortOrder] = sort([pop.Cost]);
pop = pop(sortOrder);


%% BBO Parameters    BBO�㷨����

KeepRate=0.2;                   % Keep Rate    ��������Ϣ�ر���
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats    ÿ����������Ϣ������/��������

nNew=nPop-nKeep;                % Number of New Habitats    ����Ϣ����Ŀ

% Ǩ��Ǩ��������
mu=linspace(1,0,nPop);          % Emmigration Rates    Ǩ����
%%����һ������Ϊ nPop �������������� 1 �� 0 �ȼ������ֵ��
%%mu ����Ǩ���ʣ�Emmigration rates��������ʾÿ����Ϣ�ص���Ǩ���ʣ�
%%�� BBO �㷨�У���Ǩ�������𽥼��ٵģ��������������ֵ�Ǵ� 1 �𽥼�С�� 0��
lambda=1-mu;                    % Immigration Rates    Ǩ����

%%ԭʼ��������
alpha=0.9;  %Ǩ������,������Ϣ��֮���Ǩ�Ƴ̶�
lambdaFactor = 0.9; % �������ӣ��ɸ���ʵ�����
muFactor = 1.1; % �������ӣ��ɸ���ʵ�����
pMutationInitial = 0.02;  %ͻ�����
sigma=0.02*(VarMax-VarMin); %ͻ�����
probability = 0.5;  %Ⱥ�����ƶ���ֵ


%%==========================================================================


% Best Solution Ever Found
BestSol=pop(1);

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);

%% BBO Main Loop
AvgCost = zeros(MaxIt, 1);  % ����һ�����������洢ÿһ����ƽ������

for it=1:MaxIt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%����ÿ����Ϣ�صı������
    costThreshold = median([pop.Cost]); % �ɱ���ֵ������ȡ��λ��
    for i = 1:length(pop)
        if pop(i).Cost > costThreshold
            pop(i).pMutation = pMutationInitial * 1.3; % ��Ӧ�Ȳ����Ϣ�ر����ʼ�da
        else
            pop(i).pMutation = pMutationInitial * 0.7;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newpop=pop;
    
    %%�������ƶȲ�����Ƿ���Ҫ�����������ʣ�
    similarityCount = SimilarityCalculation(pop,nPop);
    if similarityCount > nPop * 0.4  % �����Ⱥ�г��� 40% �ĸ��������
        disp('High similarity detected, applying diversity adjustment...');
        % ��������������������ԣ�������ͻ����ʻ�����Ǩ��
        pMutationInitial = 0.055;  % ����ͻ�����
        sigma = 0.03 * (VarMax - VarMin);  % ���ӱ��췶Χ

        mu = mu * 0.8;    % Ǩ����
        lambda = lambda * 1.2;
    else
        pMutationInitial = 0.02;  % �ָ�������ͻ�����
        sigma = 0.02 * (VarMax - VarMin);

    end
    EP = zeros(1, nPop);
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
            if rand<=pop(i).pMutation
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

