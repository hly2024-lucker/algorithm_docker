function Chrom = InitPop(NIND, N)
%% ��ʼ����Ⱥ����
% ���룺
%   NIND - ��Ⱥ��С�������������
%   N    - ÿ������Ļ��򳤶ȣ���Ӧ�ڳ��е�������
% �����
%   Chrom - ��ʼ������Ⱥ����ÿһ�б�ʾһ�������·��

Chrom = zeros(NIND, N);  % Ԥ������Ⱥ����
for i = 1:NIND
    Chrom(i, :) = randperm(N);  % Ϊÿ�����������������·��
end
end
