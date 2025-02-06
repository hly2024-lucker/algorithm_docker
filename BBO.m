function [Leader_score,Leader_pos,Convergence_curve1,Convergence_curve2,Convergence_curve3]=BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost1,ComputeCost2,ComputeCost3)

disp(pop(3).Cost);



%% Problem Definition  �Ż����ⶨ��

CostFunction1=@(x) ComputeCost1(x);        % Cost Function
CostFunction2=@(x) ComputeCost2(x);        % Cost Function
CostFunction3=@(x) ComputeCost3(x);        % Cost Function
% 
% Max_Cost1 = zeros(MaxIt, 1);  % ����һ�����������洢ÿһ���ĺ���1��������
% Max_Cost2 = zeros(MaxIt, 1);  % ����һ�����������洢ÿһ���ĺ���2������
% Max_Cost3 = zeros(MaxIt, 1);  % ����һ�����������洢ÿһ���ĺ���3������

Max_Cost1 = 0;  % �洢����1������
Max_Cost2 = 0;  % �洢����2������
Max_Cost3 = 0;  % �洢����3������


VarSize=[1 nVar];   % Decision Variables Matrix Size  ���߾���
VarMin = varMin(1);
VarMax = varMax(1);

%VarMin=-10;         % Decision Variables Lower Bound  ���߱������½�
%VarMax= 10;         % Decision Variables Upper Bound  ���߱������Ͻ�

%% BBO Parameters    BBO�㷨����

KeepRate=0.2;                   % Keep Rate    ��������Ϣ�ر���
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats    ÿ����������Ϣ������/��������
nNew=nPop-nKeep;                % Number of New Habitats    ����Ϣ����Ŀ

% Migration Rates
mu=linspace(1,0,nPop);          % Emmigration Rates    Ǩ����
%%����һ������Ϊ nPop �������������� 1 �� 0 �ȼ������ֵ��
%%mu ����Ǩ���ʣ�Emmigration rates��������ʾÿ����Ϣ�ص���Ǩ���ʣ�
%%�� BBO �㷨�У���Ǩ�������𽥼��ٵģ��������������ֵ�Ǵ� 1 �𽥼�С�� 0��
lambda=1-mu;                    % Immigration Rates    Ǩ����
alpha=0.9;  %Ǩ������,������Ϣ��֮���Ǩ�Ƴ̶�
pMutation=0.1;  %ͻ�����
sigma=0.02*(VarMax-VarMin); %ͻ�����
probability = 0.5;  %Ⱥ�����ƶ���ֵ

%���ֺ������ܴ��ۺ�����Ȩ��
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
%%repmat(habitat, nPop, 1) ������һ������ nPop ����Ϣ�أ�habitat �ṹ��ĸ��ƣ������飬
%%ÿ����Ϣ�ض����� Position �� Cost �����ֶ�
% pop=repmat(habitat,nPop,1);

% Initialize Habitats
for i=1:nPop
%     pop(i).Position=unifrnd(VarMin,VarMax,VarSize); %��ʼ��������Ϣ�صĻ���
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
    %��ʼ��������Ⱥ�ڸ�����ܴ���
    pop(i).Cost = m1*pop(i).Cost1 + m2*pop(i).Cost2 + m3*pop(i).Cost3;
end 
%% BBO Main Loop
AvgCost = zeros(MaxIt, 1);  % ����һ�����������洢ÿһ����ƽ������

for it=1:MaxIt
    
    newpop=pop;

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
    
    % ����ÿ��ÿ�������µ���Сֵ
    [~, SortOrder]=sort([newpop.Cost1]);
    BestCost1(it) = newpop(1).Cost1;
    [~, SortOrder]=sort([newpop.Cost2]);
    BestCost2(it) = newpop(1).Cost2;
    [~, SortOrder]=sort([newpop.Cost3]);
    BestCost3(it) = newpop(1).Cost3;
    
    for i=1:nPop
        %%���㵱ǰȺ��Ļ�ϴ���
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

    % Select Next Iteration Population  ��Ӣ����
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
    AvgCost(it) = avgCost;  % ����ÿһ�����и����ƽ������
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': bboBest Cost = ' num2str(BestCost(it)) ', Average Cost = ' num2str(avgCost)]);
    
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
Convergence_curve1 = BestCost1;
Convergence_curve2 = BestCost2;
Convergence_curve3 = BestCost3;

end

