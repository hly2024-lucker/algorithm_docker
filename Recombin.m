function SelCh = Recombin(SelCh,g,G)  % g �ǵ�ǰ������G ��������
NSel = size(SelCh,1);  % ѡ��ĸ�������
p = (1 + power(g/G,1/3)) / 3;  % ���ݵ�ǰ�������㽻�����
for i = 1:2:NSel - mod(NSel,2)
    % ����ÿ��ѡ��ĸ��壬���н�����������ڼ���ĸ���
    if is_seem(SelCh(i,:),SelCh(i+1,:)) < p
        [SelCh(i,:),SelCh(i+1,:)] = intercross(SelCh(i,:),SelCh(i+1,:));
    end
end   
