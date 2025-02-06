clear
clc
close all

citys = load('lin318.txt');  % ���س�������
citys(:,1) = [];            % ȥ����һ�У����費��Ҫ��
n = size(citys,1);          % ��������
D = zeros(n,n);             % ��ʼ���������
D = Distance(citys);        % ���ɳ���֮��ľ������

NIND = 45;  % ��Ⱥ��С������������
alpha = 1;  % ��Ϣ����Ҫ�� 
beta = 5;   % ����ʽ������Ҫ��
rho = 0.15;  % ��Ϣ�ػӷ�ϵ��
Q = 1;      % ��Ϣ�س���
Eta = 1./D; % ����ʽ��Ϣ���󣨾���ĵ�����
Tau = ones(n,n);  % ��ʼ����Ϣ�ؾ���
iter_max = 600;   % ����������
N = size(D,1);          % ���� D ��ά�ȣ����������� (���������� 34x34)
m = NIND;               % ������������Ⱥ��С��

Table = zeros(m,n);     % ·����¼��ÿֻ���ϵķ���˳��m �� n �У�


Chrom = more_init(NIND, N, D);  % ��ʼ����Ⱥ
Table = Chrom(1:NIND,:);        % �洢·���ı��
Length = PathLength(D, Chrom);  % ����·������


iter = 1;                     % ������ʼֵ
Route_best = zeros(iter_max, n);   % ÿ�ε��������·��
Length_best = zeros(iter_max, 1);  % ÿ�ε��������·������
Length_ave = zeros(iter_max, 1);   % ÿ�ε�����ƽ��·������

%%%%%%%%%%%%%%%%%% ����˫��̰�����Ƶĳ�ʼ����Ⱥ����������Ϣ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��ʼ����Ⱥ��Ϣ��
% ������ϵ�·������
    Delta_Tau = zeros(n, n);  % ��Ϣ����������
    k1=0;
    k2=0;
    min_len=min(Length);
    max_len=max(Length);
    
    %�������·�������·���ĸ���
    for i=1:m
        if min_len==Length(i)
            k1=k1+1;
        elseif max_len==Length(i)
            k2=k2+1;
        end
    end
    L_min=zeros(n,n);
    L_max=zeros(n,n);
    [min_Length,min_index]=min(Length);
    [max_Length,max_index]=max(Length);
    for i=1:(n-1)%�ֲ��������·��
        L_min(Table(min_index,i),Table(min_index,i+1))=(Q/Length(min_index))*k1;
        L_max(Table(max_index,i),Table(max_index,i+1))=(Q/Length(max_index))*k2;
    end
    L_min(Table(min_index,n),Table(min_index,1))=(Q/Length(min_index))*k1;
    L_max(Table(max_index,n),Table(max_index,1))=(Q/Length(max_index))*k2;
    %(��ǿ���·������Ϣ�أ��������·������Ϣ��)
    Tau = (1 - rho) * Tau + Delta_Tau+L_min-L_max;  % ������Ϣ�ؾ����µĹ�ʽ

%% ����Ѱ������·��
while iter <= iter_max
    % �������ÿֻ���ϵ���ʼ����
    start = zeros(m, 1);
    for i = 1:m
        temp = randperm(n);  % ����������г��б��
        start(i) = temp(1);  % ÿֻ����ѡ���һ��������Ϊ���
    end
    Table(:, 1) = start;  % ��ʼ��·����¼��ĵ�һ��Ϊ��ʼ����

    % �������г����б�
    citys_index = 1:n;

    % �������ѡ��·��
    for i = 1:m
        % Ϊÿֻ����ѡ������·��
        for j = 2:n
            tabu = Table(i, 1:(j-1));  % �ѷ��ʵĳ��У������б�
            allow_index = ~ismember(citys_index, tabu);  % ʣ��δ���ʳ���
            allow = citys_index(allow_index);  % �ɷ��ʳ��м���

            % �������δ���ʳ��еĸ���
            P = allow;
            for k = 1:length(allow)
                P(k) = Tau(tabu(end), allow(k))^alpha * Eta(tabu(end), allow(k))^beta;
            end
            P = P / sum(P);  % ��һ������

            % ʹ�����̶ķ�ѡ����һ������
            Pc = cumsum(P);  % �ۼӸ����γ�����
            target_index = find(Pc >= rand, 1);  % ���ѡ��һ������
            target = allow(target_index);  % ȷ�����ʵĳ���
            Table(i, j) = target;  % ����·����
        end
    end

    % ����ÿֻ���ϵ�·������
    Length = zeros(m, 1);
    for i = 1:m
        Route = Table(i, :);
        for j = 1:(n - 1)
            Length(i) = Length(i) + D(Route(j), Route(j + 1));
        end
        Length(i) = Length(i) + D(Route(n), Route(1));  % ���ϻص����ľ���
    end

    % ��������·�����䳤��
    if iter == 1
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min_Length;
        Length_ave(iter) = mean(Length);
        Route_best(iter, :) = Table(min_index, :);
    else
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter - 1), min_Length);
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter, :) = Table(min_index, :);
        else
            Route_best(iter, :) = Route_best(iter - 1, :);
        end
    end

    % ������Ϣ�ؾ���
    Delta_Tau = zeros(n, n);
 k1=0;
    k2=0;
    min_len=min(Length);
    max_len=max(Length);
    
    %�������·�������·���ĸ���
    for i=1:m
        if min_len==Length(i)
            k1=k1+1;
        elseif max_len==Length(i)
            k2=k2+1;
        end
    end
    L_min=zeros(n,n);
    L_max=zeros(n,n);
    [min_Length,min_index]=min(Length);
    [max_Length,max_index]=max(Length);
    for i=1:(n-1)%�ֲ��������·��
        L_min(Table(min_index,i),Table(min_index,i+1))=(Q/Length(min_index))*k1;
        L_max(Table(max_index,i),Table(max_index,i+1))=(Q/Length(max_index))*k2;
    end
    L_min(Table(min_index,n),Table(min_index,1))=(Q/Length(min_index))*k1;
    L_max(Table(max_index,n),Table(max_index,1))=(Q/Length(max_index))*k2;
    %(��ǿ���·������Ϣ�أ��������·������Ϣ��)
    Tau = (1 - rho) * Tau + Delta_Tau+L_min-L_max;  % ������Ϣ�ؾ����µĹ�ʽ


    % �����������ӣ����·����
    iter = iter + 1;
    Table = zeros(m, n);
end


%% �����ʾ
[Shortest_Length, index] = min(Length_best);  % �ҵ����·����������
Shortest_Route = Route_best(index, :);  % ��ȡ����·��

% ��ӡ����·�����䳤��
disp(['��̾���: ', num2str(Shortest_Length)]);
disp(['����·��: ', num2str([Shortest_Route Shortest_Route(1)])]);  % �ص����

%% ��ͼ
figure(1);
plot([citys(Shortest_Route, 1); citys(Shortest_Route(1), 1)], ...
     [citys(Shortest_Route, 2); citys(Shortest_Route(1), 2)], 'o-');  % ����·��
grid on;

% Ϊÿ�����б�ע���
for i = 1:size(citys, 1)
    text(citys(i, 1), citys(i, 2), ['  ' num2str(i)]);
end

% ��ע�����յ�
text(citys(Shortest_Route(1), 1), citys(Shortest_Route(1), 2), ' ���');
text(citys(Shortest_Route(end), 1), citys(Shortest_Route(end), 2), ' �յ�');

% �����������ǩ�ͱ���
xlabel('����λ�� X ��');
ylabel('����λ�� Y ��');
title(['��Ⱥ�㷨�Ż�·�� (��̾���: ' num2str(Shortest_Length) ')']);

%% ·�����ȱ仯����
figure(2);
plot(1:iter_max, Length_best, 'b', 1:iter_max, Length_ave, 'r');  % ���Ƴ��ȱ仯����
legend('��̾���', 'ƽ������');
xlabel('��������');
ylabel('����');
title('������������̾�����ƽ������仯');






