function [similarityCount] = SimilarityCalculation(pop,nPop)
% ������Ⱥ�����ƶȲ�����Ƿ񳬹���ֵ
similarityThreshold = 0.9;  % ���ƶ���ֵ
similarityCount = 0;  % ��¼���ƶȹ��ߵĸ����

for i = 1:nPop
    for j = i+1:nPop
        % ������� i �� j ֮������ƶ�
        dist = DistanceCalculation(pop(i).Position, pop(j).Position);
        % �������С��ĳ����ֵ������Ϊ��������
        if dist < similarityThreshold
            similarityCount = similarityCount + 1;
        end
    end
end


end