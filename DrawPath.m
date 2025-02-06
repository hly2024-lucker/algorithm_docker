function DrawPath(Chrom, X)
%% ����·���ĺ���
% ����:
%   Chrom - ��ʾ·���ľ���ÿһ����һ��·��
%   X     - ÿ�����е�����λ�þ��� (n x 2)��ÿ�б�ʾһ�����е� (x, y) ����

R = [Chrom(1, :) Chrom(1, 1)];  % ���·�����������ӵ�ĩβ
figure;
hold on;

% �������г��е������
plot(X(:, 1), X(:, 2), 'o', 'Color', [0.5, 0.5, 0.5]);  
% ���������
plot(X(Chrom(1, 1), 1), X(Chrom(1, 1), 2), 'rv', 'MarkerSize', 20); 

% ��עÿ�����еı��
for i = 1:size(X, 1)
    text(X(i, 1) + 0.05, X(i, 2) + 0.05, num2str(i), 'Color', [1, 0, 0]);
end

% ����·������Ӽ�ͷ
A = X(R, :);  % ��ȡ·���ϵĳ�������
for i = 2:size(A, 1)
    % ����ת���� figure ����ϵ
    [arrowx, arrowy] = dsxy2figxy(gca, A(i-1:i, 1), A(i-1:i, 2));
    % ���ƴ���ͷ��·��
    annotation('textarrow', arrowx, arrowy, 'HeadWidth', 8, 'Color', [0, 0, 1]);
end

hold off;
xlabel('������');
ylabel('������');
title('·��ͼ');
box on;

end
