% ��ʼ����Ⱥ�ĺ���
% ��ʼ����Ⱥ
function population = InitPopulation(popSize, numCities)
    population = zeros(popSize, numCities);
    for i = 1:popSize
        population(i, :) = randperm(numCities);    %����һ��(���ظ�??)������
    end
end