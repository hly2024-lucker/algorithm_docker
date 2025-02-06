function SelCh = Select(Chrom,FitnV, GGAP)
% Select 选择操作
% 输入：
% Chrom - 当前种群的个体
% FitnV - 适应度值
% GGAP - 选择比例
% 输出：
% SelCh - 选择后的个体
NIND = size(Chrom,1);
NSel = max(floor(NIND * GGAP +.5),2);
Chrlx = Sus(FitnV,NSel);
SelCh = Chrom(Chrlx,:);
end



