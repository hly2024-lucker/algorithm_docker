function SelCh = Reverse(SelCh,D)  %% D �ǳ��м�ľ������
    [row,col] = size(SelCh);  % ��ȡѡ����������������
    ObjV = PathLength(D,SelCh);  % ����·������Ӧ��ֵ
    SelCh1 = SelCh;  % ���� SelCh �ĸ���
    for i = 1:row
        r1 = randsrc(1,1,(1:col));  % ���ѡ���һ����ת��
        r2 = randsrc(1,1,(1:col));  % ���ѡ��ڶ�����ת��
        mininverse = min([r1 r2]);  % �ҵ�������ת�����Сֵ
        maxinverse = max([r1 r2]);  % �ҵ�������ת������ֵ
        SelCh1(i,mininverse:maxinverse) = SelCh1(i,maxinverse:-1:mininverse);  % ��תѡ�������ڵĳ���
    end
    ObjV1 = PathLength(D,SelCh1);  % ������·������Ӧ��ֵ
    index = ObjV1 < ObjV;  % �ҵ���Ӧ��ֵ���Ƶĸ���
    SelCh(index,:) = SelCh1(index,:);  % ������Ӧ��ֵ���õĸ���
end
