function [similarityCount] = SimilarityCalculation(pop,nPop)
% 计算种群的相似度并检查是否超过阈值
similarityThreshold = 0.9;  % 相似度阈值
similarityCount = 0;  % 记录相似度过高的个体对

for i = 1:nPop
    for j = i+1:nPop
        % 计算个体 i 和 j 之间的相似度
        dist = DistanceCalculation(pop(i).Position, pop(j).Position);
        % 如果距离小于某个阈值，则认为它们相似
        if dist < similarityThreshold
            similarityCount = similarityCount + 1;
        end
    end
end


end