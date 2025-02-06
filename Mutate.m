function SelCh = Mutate(SelCh, Distance)  % Distance 是城市之间的距离矩阵

[NSel, L] = size(SelCh);  % NSel 是选择的个体数量，L 是城市数量
NIND = size(SelCh, 1);    % NIND 是总的个体数量
D = zeros(NIND, NIND);   % 初始化距离矩阵
A = 2;                    % 常量 A
pm1 = 0.06;               % 初始突变概率
pm2 = 0.001;              % 最小突变概率
ObjV = PathLength(Distance, SelCh);  % 计算路径长度
FitnV = Fitness(ObjV);    % 计算适应度值

f_max = max(FitnV);       % 最大适应度值
f_ave = seek_ave(SelCh, FitnV);  % 计算适应度值的平均值
[TobjV, index] = sort(ObjV);  % 对路径长度进行排序并获取索引

% 计算城市之间的距离
for i = 1:NIND - 1
    for j = i + 1:NIND
        for n = 1:length(SelCh(i, :))
            D(i, j) = D(i, j) + abs(SelCh(i, n) - SelCh(j, n));
        end
        D(j, i) = D(i, j);  % 距离矩阵是对称的
    end
end

j = index(1);            % 获取适应度值最好的个体的索引
a = max(max(D));         % 距离矩阵的最大值
H = mean(D(:));          % 距离矩阵的平均值
%for i = 1:NIND - 1
    for per = 1:NSel         % 遍历每个选择的个体
        f = FitnV(per);      % 获取个体的适应度值
        if f >= f_ave        % 如果个体的适应度值大于等于平均值
            pm = pm1 - (pm1 - pm2) / (f_max - f_ave) * A * (f_max - f) * abs(D(i, j) - H) / (a - H);
        elseif f < f_ave     % 如果个体的适应度值小于平均值
            pm = pm1 * A * (D(i, j) - H) / (a - H);
        end
        if pm >= rand        % 如果随机数小于等于概率 pm，则进行突变
            R = randperm(L); % 随机生成一个城市的排列
            SelCh(per, R(1:2)) = SelCh(per, R(2:-1:1));  % 交换选中的两个城市
        end
    end
end
