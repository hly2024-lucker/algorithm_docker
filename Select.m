function SelCh = Select(Chrom,FitnV, GGAP)
% Select ѡ�����
% ���룺
% Chrom - ��ǰ��Ⱥ�ĸ���
% FitnV - ��Ӧ��ֵ
% GGAP - ѡ�����
% �����
% SelCh - ѡ���ĸ���
NIND = size(Chrom,1);
NSel = max(floor(NIND * GGAP +.5),2);
Chrlx = Sus(FitnV,NSel);
SelCh = Chrom(Chrlx,:);
end



