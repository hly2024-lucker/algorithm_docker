% % 计算成本的函数（这里是一个示例函数，需要根据具体问题来定义）
% function z=ComputeCost(x)
% 
%     z=sum(x.^2);
% 
% end


% function z = ComputeCost(x)
%     % Rastrigin 函数
%     n = numel(x); % 决策变量的数量
%     z = 10 * n + sum(x .^ 2 - 10 * cos(2 * pi * x));
% end

% function z = ComputeCost(x)
%     % Ackley 函数
%     n = numel(x);
%     term1 = -20 * exp(-0.2 * sqrt(sum(x .^ 2) / n));
%     term2 = -exp(sum(cos(2 * pi * x)) / n);
%     z = term1 + term2 + 20 + exp(1);
% end


% function z = ComputeCost(x)
%     % Griewank 函数
%     sumTerm = sum(x .^ 2) / 4000;
%     prodTerm = prod(cos(x ./ sqrt(1:numel(x))));
%     z = sumTerm - prodTerm + 1;
% end

function z = ComputeCost(x)

    % Levy 函数
    n = numel(x);
    w = 1 + (x - 1) / 4;
    term1 = sin(pi * w(1))^2;
    term2 = sum((w(1:end-1) - 1).^2 .* (1 + 10 * sin(pi * w(1:end-1) + 1).^2));
    term3 = (w(end) - 1)^2 * (1 + sin(2 * pi * w(end))^2);
    z = term1 + term2 + term3;
end

% 
% %%TSP问题
% function z = ComputeCost(x)
%     X = load('tsp225.txt');  % 加载城市坐标数据
%     X(:,1) = [];            % 移除第一列（假设是城市编号）

















