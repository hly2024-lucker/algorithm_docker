function [pop] = Remove_Calculation(pop,nPop,VarMin,VarMax,CostFunction)
% ������Ⱥ�����ƶȲ�����Ƿ񳬹���ֵ
Remove_Threshold = 0.1;  %  ȥ����ֵ

for i = 1:nPop
    for j = i+1:nPop
        % ������� i �� j ֮������ƶ�
        dist = DistanceCalculation(pop(i).Position, pop(j).Position);
        % �������С��ĳ����ֵ������Ϊ��������
        if dist < Remove_Threshold
            p = local_search(CostFunction,pop(1).Position,VarMin,VarMax);
            pcost = ComputeCost(p);
            %����ֲ�������Ľ���ţ���ȡ�����Ž�
            pop(j).Position = p;
            pop(j).Cost = pcost;
        end
        
    end
end

% % ������ƶȹ��ߣ���ִ�е�������
% if similarityCount > nPop * probability  % ���磬������� 50% �ĸ��������
%     disp('High similarity detected, applying diversity adjustment...');
%     % ����ִ��������Ӹ���ͻ���ѡ����಻ͬ�������Ǩ�ƵĲ���
% end

end