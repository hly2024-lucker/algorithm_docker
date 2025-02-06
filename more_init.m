function Chrom = InitPop(NIND, N, D)  % NIND �Ǹ���������N �ǳ�������
for i = 1:N
    City(1, i) = i;    % ���б��
end
Init = zeros(N, 4);  % ��ʼ������
for i = 1:N
    [D(i,:), id] = sort(D(i,:));  % city_id: ��¼�� i ���о�������ĳ��еı��
    for j = 1:4   % ��ȡ����� i ����� 4 �����еı��
        Init(i, j) = id(j);
    end
end  % ��ʼ�������ڽӾ������

for i = 1:NIND   % ����ÿ������
    index = zeros(1, N); % ��ʼ������ĳ��б������
    total = 1;  % ��¼��ǰ����ĳ������� 
    s = randi([1, N], 1, 1);  % ���ѡ��һ�����б�� (1~N)
    Chrom(i, 1) = s;
    index(total) = s;   % ��ѡ��ĳ��б�Ŵ�����б������
    allow = City(~ismember(City, index));  % ��ȡ��ѡ���б��
    allow(find(allow == 0)) = []; 

    while total < N
        s1 = index(total);              % ��ǰ����
        t1 = Init(s1, randi([2, 4], 1, 1)); % ���ѡ������� s1 ���ڵĳ���
        while(~ismember(t1, allow))    % ��� t1 �Ƿ��ڿ�ѡ������
            if ~ismember(Init(s1,:), allow)
                % ����ǰ����û�п�ѡ���ڽӳ���
                [m, n] = size(allow);
                a = ceil(n * rand(1, 1));
                t1 = allow(ceil(a));
            else
                % �����������ڽӳ��������ѡ��
                temp_allow = Init(s1, ismember(Init(s1,:), allow));  % ��ȡ�� s1 ��صĳ���
                [m, n] = size(temp_allow);
                a = ceil(n * rand(1, 1));
                t1 = temp_allow(a);
            end
        end
        s2 = index(1);   % ѡȡ����һ������
        t2 = Init(s2, randi([2, 4], 1, 1));
        while(~ismember(t2, allow))    % ��� t2 �Ƿ��ڿ�ѡ������
            if ~ismember(Init(s2,:), allow)  % ��ǰ���� s2 û�п�ѡ�ڽӳ���
                [m, n] = size(allow);
                a = ceil(n * rand(1, 1));
                t2 = allow(ceil(a));
            else
                % �����������ڽӳ��������ѡ��
                temp_allow = Init(s2, ismember(Init(s2,:), allow));      
                [m, n] = size(temp_allow);
                a = ceil(n * rand(1, 1));
                t2 = temp_allow(a);
            end
        end
        if D(s, s1) <= D(s, s2)    % �ж��ĸ����и���
            t = t1;
            total = total + 1;
            index(total) = t; % ���³��б��
            Chrom(i, total) = t;
        else
            t = t2;
            index(1, 2:total + 1) = index(1, 1:total);
            Chrom(i, 2:total + 1) = Chrom(i, 1:total);
            total = total + 1;
            index(1) = t;
            Chrom(i, 1) = t;
        end
        s = t;
        allow = City(~ismember(City, index));  % ���¿�ѡ����
        allow(find(allow == 0)) = []; 
    end
end
for i = 1:NIND
    p = OutputPath(Chrom(i,:)); % ���·��
end
