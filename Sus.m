function NewChrIx = Sus(FitnV,Nsel)
% SelectIndividuals 选择个体的函数
% 输入：
% FitnV - 当前种群中个体的适应度值
% Nsel - 选择的个体数量
% 输出：
% NewChrlx - 选择的个体索引
[Nind,ans]= size(FitnV);
cumfit = cumsum(FitnV);
trials = cumfit(Nind) / Nsel * (rand + (0 :Nsel - 1)');
Mf = cumfit(:, ones(1, Nsel));
Mt = trials(:, ones(1, Nind))' ;
[NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel) ; Mf(1 : Nind-1,:) ] <= Mt);
[ans, shuf] = sort(rand(Nsel, 1));
NewChrIx = NewChrIx(shuf);

end

