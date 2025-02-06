%
% ��Ȩ���� (c) 2015, Yarpiz (www.yarpiz.com)
% ��������Ȩ�������Ķ���license.txt�����˽����֤���
%
% ��Ŀ����: YPEA113
% ��Ŀ����: �Ľ��� BBO �㷨��������ǿ�ı�����Ժ�������ֵ��
% ������: Yarpiz (www.yarpiz.com)
%
% ������: S. Mostapha Kalami Heris (Yarpiz �Ŷӳ�Ա)
% �޸���: [��������]
%
function [Leader_score,Leader_pos,Convergence_curve]=local_restart_Mutate_new_BBO(pop,nPop,MaxIt,varMin,varMax,nVar,ComputeCost)
%%[Best_score,Best_pos,cg_curve]

%% ���ⶨ��
CostFunction = @(x) ComputeCost(x);    % �ɱ�����

%nVar = 5;                              % ���߱���������
VarSize = [1 nVar];                    % ���߱�������Ĵ�С
VarMin = varMin(1);                         % ���߱������½�
VarMax = varMax(1);                           % ���߱������Ͻ�

%% BBO ����
% MaxIt = 600;                          % ����������
% nPop = 50;                             % ��Ⱥ��ģ����Ϣ��������

KeepRate = 0.2;                        % ��Ӣ������
nKeep = round(KeepRate * nPop);        % ÿ����������Ϣ������
nNew = nPop - nKeep;                   % ÿ�������ɵ���Ϣ������

pmu = linspace(1, 0, nPop);             % Ǩ����
plambda = 1 - pmu;                       % Ǩ����

alpha = 0.9;                           % Ǩ������
pMutation = 0.1;                       % ��ʼ�������
sigma = 0.02 * (VarMax - VarMin);      % �������
convergenceThreshold = 1e-400;           % ������ֵ


% ���ɱ�������Ⱥ
[~, SortOrder] = sort([pop.Cost]);
pop = pop(SortOrder);

% ��ѽ������
BestSol = pop(1);
same_iter = 0;
% �洢��Ѵ��ۺ�ƽ�����۵�����
BestCost = zeros(MaxIt, 1);
AvgCost = zeros(MaxIt, 1);

%% BBO ��ѭ��

for it = 1:MaxIt
    newpop = pop;
    %% ����ӦǨ���ʵ���(Ŀǰ�Ѿ��Ƴ��˹���)
    % ���ݵ�ǰ��������Ǩ����
    pmu = linspace(1, 0, nPop)* (1 - (it / MaxIt)*0.8);    % �������СǨ����
    plambda = (1 - pmu)* (1 - (it / MaxIt)*0.8);
    mu = pmu;
    lambda = plambda;
    % ������Ⱥ���ƶȲ������������
    similarityCount = SimilarityCalculation(pop, nPop);
    if similarityCount > nPop * 0.4
        disp('��⵽�����ƶȣ���ǿ��Ⱥ������...');
        pMutation = 0.05 + 0.01 * (1 - it / (MaxIt+1));  %��ʱ����������������ƶ�+������������ڱ���ǿ��
        sigma = 0.03 * (VarMax - VarMin);  % ���ӱ��췶Χ,�������С�������
        mu = pmu*0.9;
        lambda = plambda * 1.1;
        if similarityCount > nPop * 0.6
            mu = pmu*0.8;
            lambda = plambda * 1.2;
        end
    else
        pMutation = 0.02 + 0.03 * (1 - it / (MaxIt+1));  % ������������ڱ���ǿ��
        sigma = (0.02 + 0.01 * (1 - it / (MaxIt+1)))* (VarMax - VarMin);  % �������С�������
    end

 
    %% ������������㷨ͣ��ʱ���������ѡ�񲿷ָ���  �����������������Ҳ��������Ϣ�س�ʼ�������pop(i).Position = unifrnd(VarMin, VarMax, VarSize);��
        if(mod(same_iter,MaxIt/3) == 0&&same_iter>0)
            random_restart_idx = randperm(nPop, round(nPop * 0.3));  % ���ѡ�񲿷ָ���  15��
            % ������ֵPosition
            for i = 1:length(random_restart_idx)
                t = random_restart_idx(i);
                pop(t).Position = random_restart(1,nVar,VarMin,VarMax);  % �����ֵPosition
                pop(t).Cost = ComputeCost(pop(t).Position);
            end
            same_iter = 0;
        end 



    %% Ǩ�ƺͱ���
    for i = 1:nPop

        for k = 1:nVar
            % Ǩ�Ʋ���
            if rand <= lambda(i)
            EP = mu;
            EP(i) = 0;                % �ų���ǰ��Ϣ��
            EP = EP / sum(EP);        % ��һ������
            % ѡ��Դ��Ϣ�ز�����Ǩ��
            j = RouletteWheelSelection(EP);
                newpop(i).Position(k) = pop(i).Position(k) + alpha * ( pop(j).Position(k) - pop(i).Position(k) );
            end
            % ����Ӧ���죨���и�˹������
            if rand <= pMutation
                newpop(i).Position(k) = newpop(i).Position(k) + sigma * randn;
            end
        end

        % �߽�����
        newpop(i).Position = max(newpop(i).Position, VarMin);
        newpop(i).Position = min(newpop(i).Position, VarMax);

        % ������λ�õĳɱ�
        newpop(i).Cost = CostFunction(newpop(i).Position);
    end
    
    % ������Ⱥ��������
    [~, SortOrder] = sort([newpop.Cost]);
    newpop = newpop(SortOrder);

    % ѡ����һ����Ⱥ
    pop = [pop(1:nKeep); newpop(1:nNew)];

    % ������Ⱥ��������ѽ������
    [~, SortOrder] = sort([pop.Cost]);
    pop = pop(SortOrder);

    %% �����������������Ƚϴ���ֲ��������ƽ��иĽ�,�����Ž�ֲ�̽��

    if(it>MaxIt*0.9||mod(same_iter,20)==0)
        p = local_search(CostFunction,pop(1).Position,VarMin,VarMax);
        pcost = ComputeCost(p);
        %����ֲ�������Ľ���ţ���ȡ�����Ž�
        if(pcost<pop(1).Cost)
            pop(1).Position = p;
            pop(1).Cost = pcost;
        end
    end
    BestSol = pop(1);
    % �洢��Ѻ�ƽ������
    BestCost(it) = BestSol.Cost;
    AvgCost(it) = mean([pop.Cost]);
    
    %%%������䣬���е���
    tmp = max(it-1,1);
%     if it > 1 && abs(BestCost(tmp) - BestSol.Cost) < abs(BestCost(tmp))*(1/MaxIt)
    if(abs(BestCost(tmp) - BestSol.Cost)==0)
        disp('---------------------------------');
        same_iter = same_iter+1;
        if(same_iter>50)
            pMutation = 0.05 + 0.01 * (1 - it / (MaxIt+1));
            sigma = 0.03 * (VarMax - VarMin);  % ���ӱ��췶Χ,�������С�������
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

    % ��ʾ������Ϣ
    disp(['���� ' num2str(it) ': ��Ѵ��� = ' num2str(BestCost(it)) ', ƽ������ = ' num2str(AvgCost(it))]);

    % �������
    if it > 1 && abs(BestCost(it) - BestCost(it - 1)) < convergenceThreshold
        disp('�ﵽ������ֵ����ǰֹͣ...');
        BestCost = BestCost(1:it);
        AvgCost = AvgCost(1:it);
        break;
    end
end

% %% ���
% figure;
% semilogy(BestCost, 'LineWidth', 2);
% hold on;
% semilogy(AvgCost, 'LineWidth', 2, 'LineStyle', '-');
% xlabel('��������');
% ylabel('����');
% grid on;
% legend('��Ѵ���', 'ƽ������');

Leader_pos =  BestSol.Position;
Leader_score =  BestSol.Cost;
Convergence_curve = BestCost;
   %%[Best_score,Best_pos,cg_curve]  
%[Leader_pos,Leader_score,Convergence_curve] = [BestSol.Cost,BestSol.Position,BestCost];
end
