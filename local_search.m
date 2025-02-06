% 局部搜索操作 - 梯度下降式优化
%obj_func 目标函数（求解最小值）
%solution 当前最优解、待优化的解
%lb 变量最小值，ub 最大值
%local_search_iters 局部搜索次数
function best_solution = local_search(obj_func,solution,lb,ub)
    local_search_iters = 50;
    alpha = 0.01;  % 学习率，控制更新步长
    best_solution = solution;

    best_fitness = obj_func(solution);

    for i = 1:local_search_iters
        % 计算当前解的梯度（近似）
        grad = numerical_gradient(obj_func, best_solution);
        
        % 梯度下降更新：沿梯度方向反向更新解
        new_solution = best_solution - alpha * grad;
        
        % 修正位置，确保解在边界内
        new_solution = max(min(new_solution, ub), lb);
        
        % 计算新解的适应度
        new_fitness = obj_func(new_solution);
        
        % 如果新解更好，则更新解
        if new_fitness < best_fitness
            best_solution = new_solution;
            best_fitness = new_fitness;
        end
    end
end