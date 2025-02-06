clc
clear
close all
%%
nPop=30; % 种群数

Max_iter = 1000; % 最大迭代次数

dim = 50; % 可选 2, 10, 30, 50, 100

%%  选择函数

Function_name1=1 ; % 函数名： 1 - 30
[lb1,ub1,dim,fobj1] = Get_Functions_cec2017(Function_name1,dim);

Function_name2= 12; % 函数名： 1 - 30  (12,17,19)
[lb2,ub2,dim,fobj2] = Get_Functions_cec2017(Function_name2,dim);

Function_name3=21; % 函数名： 1 - 30
[lb3,ub3,dim,fobj3] = Get_Functions_cec2017(Function_name3,dim);

lb = max(max(lb1,lb2),lb3);
ub = min(min(ub1,ub2),ub3);

nVar = dim;
VarSize=[1 nVar];   % Decision Variables Matrix Size  决策矩阵 
VarMin = lb(1);
VarMax = ub(1);

%% BBO Parameters    BBO算法参数

MaxIt=Max_iter;          % Maximum Number of Iterations 最大迭代次数

%% Initializationxl

% Empty Habitat
habitat.Position=[];
habitat.Cost1=[];
habitat.Cost2=[];
habitat.Cost3=[];
habitat.Cost=[];

% Create Habitats Array
%%repmat(habitat, nPop, 1) 会生成一个包含 nPop 个栖息地（habitat 结构体的复制）的数组，
%%每个栖息地都包含 Position 和 Cost 两个字段
pop=repmat(habitat,nPop,1);

% Initialize Habitats初始化各栖息地的代价
for i=1:nPop
%     pop(i).Position=unifrnd(VarMin,VarMax,VarSize); %初始化各个栖息地的环境
%     pop(i).Cost=fobj(pop(i).Position);
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize); %初始化各个栖息地的环境
    
    pop(i).Cost1=fobj1(pop(i).Position);
    pop(i).Cost2=fobj2(pop(i).Position);
    pop(i).Cost3=fobj3(pop(i).Position);
end


%% 调用算法
tic

[Best_score1,Best_pos1,cg_curve11,cg_curve12,cg_curve13]=BBO(pop,nPop,Max_iter,lb,ub,dim,fobj1,fobj2,fobj3);
%[Best_score2,Best_pos2,cg_curve21,cg_curve22,cg_curve23]=change3_mu_lambda_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj1,fobj2,fobj3);
% [Best_score3,Best_pos3,cg_curve3]=mu_lambda_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score4,Best_pos4,cg_curve4] = change3_mu_lambda_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
[Best_score5,Best_pos5,cg_curve51,cg_curve52,cg_curve53] = change4_mu_lambda_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj1,fobj2,fobj3);
%[Best_score6,Best_pos6,cg_curve6] = change5_mu_lambda_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj1,fobj2,fobj3);



% [Best_score2,Best_pos2,cg_curve2]=new_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score3,Best_pos3,cg_curve3]=Mutate_new_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% % [Best_score4,Best_pos4,cg_curve4] = restart_Mutate_new_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score5,Best_pos5,cg_curve5] = local_restart_Mutate_new_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score6,Best_pos6,cg_curve6] = change_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score7,Best_pos7,cg_curve7] = new_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score8,Best_pos8,cg_curve8] = change3_orignal_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);
% [Best_score9,Best_pos9,cg_curve9] = change2_BBO(pop,nPop,Max_iter,lb,ub,dim,fobj);


toc

%% plot
figure('Position',[400 200 300 250])
semilogy(cg_curve11,'Color','r','Linewidth',2);
hold on;  % 继续在同一图上绘制
% semilogy(cg_curve21,'Color','b','Linewidth',2);
% hold on;  % 继续在同一图上绘制
% semilogy(cg_curve3,'Color','g','Linewidth',2);
% hold on;  % 继续在同一图上绘制
% semilogy(cg_curve4,'Color','y','Linewidth',2);
% hold on;  % 继续在同一图上绘制
semilogy(cg_curve51,'Color','k','Linewidth',2);

title(['Convergence curve1, Dim=' num2str(dim)])
xlabel('Iteration');
% ylabel(['Best score F' num2str(Function_name1) ]);
axis tight
grid on
box on
set(gca,'color','none')
legend('BBO','IBBO');

% 
% 
figure('Position',[600 400 500 450])
semilogy(cg_curve12,'Color','r','Linewidth',2);
hold on;  % 继续在同一图上绘制
% semilogy(cg_curve22,'Color','b','Linewidth',2);
% hold on;  % 继续在同一图上绘制
% semilogy(cg_curve3,'Color','g','Linewidth',2);
% hold on;  % 继续在同一图上绘制
% semilogy(cg_curve4,'Color','y','Linewidth',2);
% hold on;  % 继续在同一图上绘制
semilogy(cg_curve52,'Color','k','Linewidth',2);


title(['Convergence curve2, Dim=' num2str(dim)])
xlabel('Iteration');
% ylabel(['Best score F' num2str(Function_name1) ]);
axis tight
grid on
box on
set(gca,'color','none')
legend('BBO','IBBO');


figure('Position',[800 300 400 350])
semilogy(cg_curve13,'Color','r','Linewidth',2);
hold on;  % 继续在同一图上绘制
% semilogy(cg_curve23,'Color','b','Linewidth',2);
% hold on;  % 继续在同一图上绘制
% semilogy(cg_curve3,'Color','g','Linewidth',2);
% hold on;  % 继续在同一图上绘制
% semilogy(cg_curve4,'Color','y','Linewidth',2);
% hold on;  % 继续在同一图上绘制
semilogy(cg_curve53,'Color','k','Linewidth',2);



title(['Convergence curve3, Dim=' num2str(dim)])
xlabel('Iteration');
% ylabel(['Best score F' num2str(Function_name1) ]);
axis tight
grid on
box on
set(gca,'color','none')
legend('BBO','IBBO');





% legend('BBO', 'Similar_BBO','Mutate_BBO', 'Restart_BBO','LocalSearch_BBO','change_BBO');
% legend('BBO', 'mu_lambda_BBO','change3_mu_lambda_BBO','change4_mu_lambda_BBO');
% legend('BBO','new_BBO', 'mu_lambda_BBO','change3_mu_lambda_BBO','change4_mu_lambda_BBO');
