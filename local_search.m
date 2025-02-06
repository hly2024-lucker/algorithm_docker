% �ֲ��������� - �ݶ��½�ʽ�Ż�
%obj_func Ŀ�꺯���������Сֵ��
%solution ��ǰ���Ž⡢���Ż��Ľ�
%lb ������Сֵ��ub ���ֵ
%local_search_iters �ֲ���������
function best_solution = local_search(obj_func,solution,lb,ub)
    local_search_iters = 50;
    alpha = 0.01;  % ѧϰ�ʣ����Ƹ��²���
    best_solution = solution;

    best_fitness = obj_func(solution);

    for i = 1:local_search_iters
        % ���㵱ǰ����ݶȣ����ƣ�
        grad = numerical_gradient(obj_func, best_solution);
        
        % �ݶ��½����£����ݶȷ�������½�
        new_solution = best_solution - alpha * grad;
        
        % ����λ�ã�ȷ�����ڱ߽���
        new_solution = max(min(new_solution, ub), lb);
        
        % �����½����Ӧ��
        new_fitness = obj_func(new_solution);
        
        % ����½���ã�����½�
        if new_fitness < best_fitness
            best_solution = new_solution;
            best_fitness = new_fitness;
        end
    end
end