function SelCh = Recombin(SelCh,g,G)  % g 是当前代数，G 是最大代数
NSel = size(SelCh,1);  % 选择的个体数量
p = (1 + power(g/G,1/3)) / 3;  % 根据当前代数计算交叉概率
for i = 1:2:NSel - mod(NSel,2)
    % 对于每对选择的个体，进行交叉操作，基于计算的概率
    if is_seem(SelCh(i,:),SelCh(i+1,:)) < p
        [SelCh(i,:),SelCh(i+1,:)] = intercross(SelCh(i,:),SelCh(i+1,:));
    end
end   
