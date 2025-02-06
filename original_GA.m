clear
clc
close all

X = load('bier127.txt');  % ���س�����������
X(:,1) = [];            % �Ƴ���һ�У������ǳ��б�ţ�

% 33958.5975 //0.02

NIND = 45;             % ��Ⱥ��ģ
MAXGEN = 8000;          % ������
Pc = 0.9;               % �������
Pm = 0.02;              % �������
GGAP = 0.9;             % ѡ�����
D = Distance(X);        % ���ɾ������
N = size(D,1);          % ��ȡ��������������Ϊ34��

%% ��ʼ����Ⱥ
Chrom = InitPop(NIND, N);

%% �������·��ͼ
DrawPath(Chrom(1,:), X)
pause(0.0001);

%% ������·�������ܾ���
disp('��ʼ����Ⱥ�е�һ�����·��:');
OutputPath(Chrom(1,:));
Rlength = PathLength(D, Chrom(1,:));
disp(['�ܾ���:', num2str(Rlength)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%% �Ż�����
gen = 0;
figure;
hold on; box on;
xlim([0, MAXGEN]);
title('�Ż�����');
xlabel('����');
ylabel('����·������');
ObjV = PathLength(D, Chrom);  % ����·������
preObjV = min(ObjV);
max_N = MAXGEN / 8;  % ����������ͽ����ʵ�������
total = 1;
tmp = 0;
while gen < MAXGEN
    %% ������Ӧ��
    ObjV = PathLength(D, Chrom);  % ����·������
    line([gen-1, gen], [preObjV, min(ObjV)]);
    pause(0.0001);
    preObjV = min(ObjV);
    % ���������ϴ������Ӧ��
    preFitnV = tmp; 
    FitnV = Fitness(ObjV);
    tmp = max(FitnV);
    if tmp <= preFitnV
        total = total + 1;
        fprintf('+1\n');
        if total == max_N
            fprintf("����ͣ�ͣ�ִ����һ�㷨\n");
            break;
        end
    else
        total = 1;
    end

    %% ѡ�����
    SelCh = Select(Chrom, FitnV, GGAP);
    % �������
    SelCh = original_Recombin(SelCh, Pc);
    %% �������
    SelCh = original_Mutate(SelCh, Pm);
    % %% ��ת���� (��ѡ)
    % SelCh = Reverse(SelCh, D);
    %% �ز����Ӵ�������Ⱥ
    Chrom = Reins(Chrom, SelCh, ObjV);
    %% ���´���
    gen = gen + 1;
end

%% ��������·��ͼ
ObjV = PathLength(D, Chrom);  % ����·������
[minObjV, minInd] = min(ObjV);
Shortest_Route = Chrom(minInd(1), :);

%% ��ͼ
figure(3);
plot([X(Shortest_Route, 1); X(Shortest_Route(1), 1)], ...
    [X(Shortest_Route, 2); X(Shortest_Route(1), 2)], 'o-');
grid on;
for i = 1:size(X, 1)
    text(X(i, 1), X(i, 2), ['  ' num2str(i)]);
end
text(X(Shortest_Route(1), 1), X(Shortest_Route(1), 2), '      ���');
text(X(Shortest_Route(end), 1), X(Shortest_Route(end), 2), '      �յ�');
xlabel('����λ�ú�����');
ylabel('����λ��������');
title(['���ڴ�ͳ�Ŵ��㷨��·���Ż�(��̾��룺' num2str(minObjV) ')']);

%% �������·�����ܾ���
disp('����·����');
p = OutputPath(Chrom(minInd(1), :));
disp(['�ܾ���:', num2str(ObjV(minInd(1)))]); 
disp('-----------------------------------------------');
