% 初始化种群函数
function population = random_restart(population_size, dimensions,VarMin,VarMax)
    % 生成一个初始化种群，个体在 [-5, 5] 范围内随机生成
    %生成[0,1]矩阵然后乘十，减去五可以得到【-5，5】
    scope = VarMax - VarMin;
    population = (rand(population_size, dimensions) * scope) - abs(VarMin);
end







