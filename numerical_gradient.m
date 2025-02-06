% ��ֵ�ݶȼ��㣨���޲�ַ���
%obj_funcĿ�꺯��
%solution ��ǰ��
function grad = numerical_gradient(obj_func, solution)
    epsilon = 1e-6;  % ΢С�Ŷ�
    grad = zeros(size(solution));
    for i = 1:length(solution)
        perturbed_solution = solution;
        perturbed_solution(i) = perturbed_solution(i) + epsilon;
        grad(i) = (obj_func(perturbed_solution) - obj_func(solution)) / epsilon;
    end
end