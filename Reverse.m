function SelCh = Reverse(SelCh,D)  %% D 是城市间的距离矩阵
    [row,col] = size(SelCh);  % 获取选择个体的行数和列数
    ObjV = PathLength(D,SelCh);  % 计算路径的适应度值
    SelCh1 = SelCh;  % 创建 SelCh 的副本
    for i = 1:row
        r1 = randsrc(1,1,(1:col));  % 随机选择第一个反转点
        r2 = randsrc(1,1,(1:col));  % 随机选择第二个反转点
        mininverse = min([r1 r2]);  % 找到两个反转点的最小值
        maxinverse = max([r1 r2]);  % 找到两个反转点的最大值
        SelCh1(i,mininverse:maxinverse) = SelCh1(i,maxinverse:-1:mininverse);  % 反转选定区间内的城市
    end
    ObjV1 = PathLength(D,SelCh1);  % 计算新路径的适应度值
    index = ObjV1 < ObjV;  % 找到适应度值改善的个体
    SelCh(index,:) = SelCh1(index,:);  % 更新适应度值更好的个体
end
