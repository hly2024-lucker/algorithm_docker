% 初始化种群的函数
% 初始化种群
function population = InitPopulation(popSize, numCities)
    population = zeros(popSize, numCities);
    for i = 1:popSize
        population(i, :) = randperm(numCities);    %返回一个(不重复??)的排列
    end
end