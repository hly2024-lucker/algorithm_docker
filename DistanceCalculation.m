% 计算两个个体的欧几里得距离
function dist = DistanceCalculation(x1, x2)
    dist = sqrt(sum((x1 - x2).^2));
end
