function Chrom = Reins(Chrom,SelCh,ObjV)
%% 更新个体种群
% 输入：
% Chrom      当前种群的个体
% SelCh      选择的个体
% ObjV       当前种群的适应度值
% 输出：
% Chrom      更新后的个体种群


NIND = size(Chrom, 1);  % 当前种群个体的数量
NSel = size(SelCh, 1);  % 被选择的个体数量
[TobjV, index] = sort(ObjV);    % 对适应度值进行排序，并返回排序后的值和索引
%% 选择适应度值最好的个体，并将其加入新的种群
Chrom = [Chrom(index(1:NIND-NSel), :); SelCh];  % 更新种群，保留适应度值最好的个体

end

