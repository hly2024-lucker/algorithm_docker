function Chrom = Reins(Chrom,SelCh,ObjV)
%% ���¸�����Ⱥ
% ���룺
% Chrom      ��ǰ��Ⱥ�ĸ���
% SelCh      ѡ��ĸ���
% ObjV       ��ǰ��Ⱥ����Ӧ��ֵ
% �����
% Chrom      ���º�ĸ�����Ⱥ


NIND = size(Chrom, 1);  % ��ǰ��Ⱥ���������
NSel = size(SelCh, 1);  % ��ѡ��ĸ�������
[TobjV, index] = sort(ObjV);    % ����Ӧ��ֵ�������򣬲�����������ֵ������
%% ѡ����Ӧ��ֵ��õĸ��壬����������µ���Ⱥ
Chrom = [Chrom(index(1:NIND-NSel), :); SelCh];  % ������Ⱥ��������Ӧ��ֵ��õĸ���

end

