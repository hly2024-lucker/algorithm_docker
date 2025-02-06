% 数值梯度计算（有限差分法）
%obj_func目标函数
%solution 当前解
function grad = numerical_gradient(obj_func, solution)
    epsilon = 1e-6;  % 微小扰动
    grad = zeros(size(solution));
    for i = 1:length(solution)
        perturbed_solution = solution;
        perturbed_solution(i) = perturbed_solution(i) + epsilon;
        grad(i) = (obj_func(perturbed_solution) - obj_func(solution)) / epsilon;
    end
end