function NewChrIx = Sus(FitnV,Nsel)
% SelectIndividuals ѡ�����ĺ���
% ���룺
% FitnV - ��ǰ��Ⱥ�и������Ӧ��ֵ
% Nsel - ѡ��ĸ�������
% �����
% NewChrlx - ѡ��ĸ�������
[Nind,ans]= size(FitnV);
cumfit = cumsum(FitnV);
trials = cumfit(Nind) / Nsel * (rand + (0 :Nsel - 1)');
Mf = cumfit(:, ones(1, Nsel));
Mt = trials(:, ones(1, Nind))' ;
[NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel) ; Mf(1 : Nind-1,:) ] <= Mt);
[ans, shuf] = sort(rand(Nsel, 1));
NewChrIx = NewChrIx(shuf);

end

