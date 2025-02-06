function Chrom = InitPop(NIND, N, D)  % NIND 是个体数量，N 是城市数量
for i = 1:N
    City(1, i) = i;    % 城市编号
end
Init = zeros(N, 4);  % 初始化数组
for i = 1:N
    [D(i,:), id] = sort(D(i,:));  % city_id: 记录与 i 城市距离最近的城市的编号
    for j = 1:4   % 获取与城市 i 最近的 4 个城市的编号
        Init(i, j) = id(j);
    end
end  % 初始化城市邻接矩阵完成

for i = 1:NIND   % 遍历每个个体
    index = zeros(1, N); % 初始化个体的城市编号数组
    total = 1;  % 记录当前个体的城市数量 
    s = randi([1, N], 1, 1);  % 随机选择一个城市编号 (1~N)
    Chrom(i, 1) = s;
    index(total) = s;   % 将选择的城市编号存入城市编号数组
    allow = City(~ismember(City, index));  % 获取可选城市编号
    allow(find(allow == 0)) = []; 

    while total < N
        s1 = index(total);              % 当前城市
        t1 = Init(s1, randi([2, 4], 1, 1)); % 随机选择与城市 s1 相邻的城市
        while(~ismember(t1, allow))    % 检查 t1 是否在可选城市中
            if ~ismember(Init(s1,:), allow)
                % 若当前城市没有可选的邻接城市
                [m, n] = size(allow);
                a = ceil(n * rand(1, 1));
                t1 = allow(ceil(a));
            else
                % 否则从允许的邻接城市中随机选择
                temp_allow = Init(s1, ismember(Init(s1,:), allow));  % 获取与 s1 相关的城市
                [m, n] = size(temp_allow);
                a = ceil(n * rand(1, 1));
                t1 = temp_allow(a);
            end
        end
        s2 = index(1);   % 选取的另一个城市
        t2 = Init(s2, randi([2, 4], 1, 1));
        while(~ismember(t2, allow))    % 检查 t2 是否在可选城市中
            if ~ismember(Init(s2,:), allow)  % 当前城市 s2 没有可选邻接城市
                [m, n] = size(allow);
                a = ceil(n * rand(1, 1));
                t2 = allow(ceil(a));
            else
                % 否则从允许的邻接城市中随机选择
                temp_allow = Init(s2, ismember(Init(s2,:), allow));      
                [m, n] = size(temp_allow);
                a = ceil(n * rand(1, 1));
                t2 = temp_allow(a);
            end
        end
        if D(s, s1) <= D(s, s2)    % 判断哪个城市更近
            t = t1;
            total = total + 1;
            index(total) = t; % 更新城市编号
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
        allow = City(~ismember(City, index));  % 更新可选城市
        allow(find(allow == 0)) = []; 
    end
end
for i = 1:NIND
    p = OutputPath(Chrom(i,:)); % 输出路径
end
