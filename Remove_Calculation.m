function [pop] = Remove_Calculation(pop,nPop,VarMin,VarMax,CostFunction)
% 计算种群的相似度并检查是否超过阈值
Remove_Threshold = 0.1;  %  去重阈值

for i = 1:nPop
    for j = i+1:nPop
        % 计算个体 i 和 j 之间的相似度
        dist = DistanceCalculation(pop(i).Position, pop(j).Position);
        % 如果距离小于某个阈值，则认为它们相似
        if dist < Remove_Threshold
            p = local_search(CostFunction,pop(1).Position,VarMin,VarMax);
            pcost = ComputeCost(p);
            %如果局部搜索后的解更优，则取代最优解
            pop(j).Position = p;
            pop(j).Cost = pcost;
        end
        
    end
end

% % 如果相似度过高，则执行调整策略
% if similarityCount > nPop * probability  % 例如，如果超过 50% 的个体对相似
%     disp('High similarity detected, applying diversity adjustment...');
%     % 可以执行类似添加更多突变或选择更多不同个体进行迁移的策略
% end

end