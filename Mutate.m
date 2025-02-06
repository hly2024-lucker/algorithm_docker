function SelCh = Mutate(SelCh, Distance)  % Distance �ǳ���֮��ľ������

[NSel, L] = size(SelCh);  % NSel ��ѡ��ĸ���������L �ǳ�������
NIND = size(SelCh, 1);    % NIND ���ܵĸ�������
D = zeros(NIND, NIND);   % ��ʼ���������
A = 2;                    % ���� A
pm1 = 0.06;               % ��ʼͻ�����
pm2 = 0.001;              % ��Сͻ�����
ObjV = PathLength(Distance, SelCh);  % ����·������
FitnV = Fitness(ObjV);    % ������Ӧ��ֵ

f_max = max(FitnV);       % �����Ӧ��ֵ
f_ave = seek_ave(SelCh, FitnV);  % ������Ӧ��ֵ��ƽ��ֵ
[TobjV, index] = sort(ObjV);  % ��·�����Ƚ������򲢻�ȡ����

% �������֮��ľ���
for i = 1:NIND - 1
    for j = i + 1:NIND
        for n = 1:length(SelCh(i, :))
            D(i, j) = D(i, j) + abs(SelCh(i, n) - SelCh(j, n));
        end
        D(j, i) = D(i, j);  % ��������ǶԳƵ�
    end
end

j = index(1);            % ��ȡ��Ӧ��ֵ��õĸ��������
a = max(max(D));         % �����������ֵ
H = mean(D(:));          % ��������ƽ��ֵ
%for i = 1:NIND - 1
    for per = 1:NSel         % ����ÿ��ѡ��ĸ���
        f = FitnV(per);      % ��ȡ�������Ӧ��ֵ
        if f >= f_ave        % ����������Ӧ��ֵ���ڵ���ƽ��ֵ
            pm = pm1 - (pm1 - pm2) / (f_max - f_ave) * A * (f_max - f) * abs(D(i, j) - H) / (a - H);
        elseif f < f_ave     % ����������Ӧ��ֵС��ƽ��ֵ
            pm = pm1 * A * (D(i, j) - H) / (a - H);
        end
        if pm >= rand        % ��������С�ڵ��ڸ��� pm�������ͻ��
            R = randperm(L); % �������һ�����е�����
            SelCh(per, R(1:2)) = SelCh(per, R(2:-1:1));  % ����ѡ�е���������
        end
    end
end
