function SelCh = Recombin(SelCh,Pc)
%%����
%Selch  ��ѡ�����
%Pc    �������
%SelCh  ��������

NSel = size(SelCh,1);
for i = 1:2:NSel - mod(NSel,2)
    if Pc>=rand   %�������
        [SelCh(i,:),SelCh(i+1,:)] = intercross(SelCh(i,:),SelCh(i+1,:));
    end
end