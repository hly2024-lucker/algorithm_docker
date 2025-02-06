clear
clc
close all

citys = load('bier127.txt');  % ���س�����������
citys(:,1) = [];            % �Ƴ���һ�У������ǳ��б�ţ�
n = size(citys,1);          % ��������
D = zeros(n,n);             
NIND = 45;                  % ��Ⱥ��С
D = Distance(citys);        % ���ɾ������
N = size(D,1);              % ��ȡ��������������Ϊ34��

m = 45;                     % ��������
alpha = 1;                  % ��Ϣ����Ҫ������
beta = 5;                   % ����������Ҫ��Ȩ��
rho = 0.15;                  % ��Ϣ�ػӷ�����
Q = 1;                      % ��ϵ��
Eta = 1./D;                 % ��������
Tau = ones(n,n);            % ��Ϣ�ؾ���
Table = zeros(m,n);         % ·����¼��
iter = 1;                   % ������ʼֵ
iter_max = 400;             % ����������
Route_best = zeros(iter_max,n);     % ÿ������·��
Length_best = zeros(iter_max,1);    % ÿ������·������
Length_ave = zeros(iter_max,1);     % ÿ��·����ƽ������


%% ������������·��
while iter <= iter_max
    % �������ÿ�����ϵ�������
    start = zeros(m,1);
    for i = 1:m
        temp = randperm(n);
        start(i) = temp(1);
    end
    Table(:,1) = start;
    
    % ����һ���հ׿ռ�
    citys_index = 1:n;
    
    % ÿֻ����·��ѡ��
    for i = 1:m
        % ÿ������·��ѡ��
        for j = 2:n
            tabu = Table(i,1:(j-1));    % �ѷ��ʵĳ��м��ϣ����ɱ�
            allow_index = ~ismember(citys_index,tabu);
            allow = citys_index(allow_index);   % �ɷ��ʵĳ��м���
            P = allow;
            
            % ������м�ת�Ƹ���
            for k = 1:length(allow)
                P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
            end
            P = P / sum(P);
            
            % ���̶ķ�ѡ����һ�����ʳ���
            Pc = cumsum(P); % �ۼ�
            target_index = find(Pc >= rand);
            target = allow(target_index(1));
            Table(i,j) = target;
        end
    end

    Length = zeros(m,1);
    for i = 1:m
        Route = Table(i,:);
        for j = 1:(n-1)
            Length(i) = Length(i) + D(Route(j), Route(j+1));
        end
        Length(i) = Length(i) + D(Route(n), Route(1));
    end

    % �������·�����뼰ƽ������
    if iter == 1
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min_Length;
        Length_ave(iter) = mean(Length);
        Route_best(iter,:) = Table(min_index,:);
    else
        [min_Length, min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter-1), min_Length);
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter,:) = Table(min_index,:);
        else
            Route_best(iter,:) = Route_best(iter-1,:);
        end
    end

    % ������Ϣ��
    Delta_Tau = zeros(n,n);
    for i = 1:m
        for j = 1:(n-1)
            Delta_Tau(Table(i,j), Table(i,j+1)) = Delta_Tau(Table(i,j), Table(i,j+1)) + Q / Length(i);
        end
        Delta_Tau(Table(i,n), Table(i,1)) = Delta_Tau(Table(i,n), Table(i,1)) + Q / Length(i);
    end
    Tau = (1 - rho) * Tau + Delta_Tau;

    % ����������һ�����·����¼��
    iter = iter + 1;
    Table = zeros(m, n);
end

%% �����ʾ
[Shortest_Length, index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['��̾��룺', num2str(Shortest_Length)]);
disp(['���·����', num2str([Shortest_Route Shortest_Route(1)])]);

%% ��ͼ
figure(1)
plot([citys(Shortest_Route,1); citys(Shortest_Route(1),1)], ...
    [citys(Shortest_Route,2); citys(Shortest_Route(1),2)], 'o-');
grid on
for i = 1:size(citys, 1)
    text(citys(i, 1), citys(i, 2), ['  ' num2str(i)]);
end
text(citys(Shortest_Route(1),1), citys(Shortest_Route(1),2), '      ���');
text(citys(Shortest_Route(end),1), citys(Shortest_Route(end),2), '      �յ�');
xlabel('����λ�ú�����');
ylabel('����λ��������');
title(['��Ⱥ�㷨·��(��̾��룺' num2str(Shortest_Length) ')']);

figure(2)
plot(1:iter_max, Length_best, 'b', 1:iter_max, Length_ave, 'r');
legend('��̾���', 'ƽ������');
xlabel('��������');
ylabel('����');
title('������̾�����ƽ������Ա�');
