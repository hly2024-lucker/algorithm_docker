Map1=[82,7;91,38;83,46;71,44;64,60;68,58;
83,69;87,76;74,78;71,71];%具体城市坐标
MaxIter=200;   %最大迭代次数
Pop_Size=150;  %种群规模
pc=0.8;        %选择概率
pm=0.1;        %变异概率
[opt,fval]=tspga(Map1,MaxIter,Pop_Size,pm,pc);

function f=fitness(fmin,fmax,froad)
f=1-(froad-fmin)/(fmax-fmin);       %适应度函数
end

function d=distance(A,B)
d=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2); %距离公式
end

function  elim=eliminate(x,y)
for n=1:length(y)
    x=x(find(x~=y(n))); %去除x中的y
end
elim=x;
end


%遗传算法计算TSP
function  [opt,fval]=tspga(Map,MaxIter,SizeScale,pm,pc)
n=max(size(Map));
DistMatrix=zeros(n,n);%初始化距离矩阵DistMatrix()，用于存储两城市之间距离
for i=1:n
    for j=1:n
        DistMatrix(i,j)=distance(Map(i,:),Map(j,:));%计算两城市之间距离，存储于DistMatrix()矩阵中
    end
end

%%%%%%%%%生成初始种群%%%%%%%%%%%
Road=ones(SizeScale,n);%初始化路径矩阵Road()
for i=1:SizeScale
    Road(i,:)=randperm(n);%随机生成初始种群（路径矩阵）
end
iter=1;
MinestRoad_fval=ones(MaxIter,1);%初始化最短里程矩阵MinestRoad_fval()的历史记录值
MinestRoad_opt=ones(MaxIter,n);%初始化最短里程路径矩阵MinestRoad_opt()的历史记录值
while iter<=MaxIter  
    Dist=zeros(SizeScale,1);%初始化里程矩阵Dist，用于存储每条路径的里程值
%%%%%%计算每条路径的里程%%%%%
    for i=1:SizeScale
    for j=1:(n-1)
       Dist(i)=Dist(i)+DistMatrix(Road(i,j),Road(i,j+1));
    end  
       Dist(i)=Dist(i)+DistMatrix(Road(i,1),Road(i,n));%计算每条路径的里程存储于Dist()矩阵中
end 
%%%%%%%%计算每条路径的适应度值%%%%%%%%%
fitmatrix=ones(SizeScale,1);
[MinRoad A]=min(Dist(:,1));%计算出最小里程值
MaxRoad=max(Dist(:,1)); %计算出最大里程值
for i=1:SizeScale
     fitmatrix(i)=fitness(MinRoad,MaxRoad,Dist(i));%计算每条路径的适应度的值，存储于fitmatrix()中
end
%%%%%%%%%%选择操作%%%%%%%%%             
[c p]=sort(fitmatrix(:,1));%对适应度值进行升序排列，c中存放升序排列适应度的值，
                           %p中存放的是c对fitmatrix()的索引，即适应度值对应的fitmatrix()中的第几行
change=20;%选出适应度值最小路径数目
for i=1:change
    Road(p(i),:)=Road(p(SizeScale),:);%选出适应度值最小的20条路径，用适应度值最大的路径替换
end
Roadnew=Road;
%%%%%%%%%%%%交叉操作%%%%%%%%%%%%
for i=1:SizeScale
u=randi([1 SizeScale],2,1);
s=u(1);
t=u(2);
if rand(1)<pc %判断是否进行交叉操作 
oldp1=Road(s,:);%随机选取两个父代染色体
oldp2=Road(t,:);
u=randi([1 n],2,1);
crossj1=u(1);%随机选取两个切点
crossj2=u(2);
minjcross=min(crossj1,crossj2);
maxjcross=max(crossj1,crossj2);
segment1=oldp1(minjcross:maxjcross);%选中的路径（染色体）片段，文中的AXY
segment2=oldp2(minjcross:maxjcross);%选中的路径（染色体）片段，文中的BXY 
oldp12=eliminate(oldp1,segment2);%在oldp1删除与segment2相同的元素
oldp21=eliminate(oldp2,segment1);%在oldp2删除与segment1相同的元素
newp1=[segment2,oldp12];%生成新路径（子代染色体）
newp2=[segment1,oldp21];%生成新路径（子代染色体）
Roadnew(s,:)=newp1;%更新种群
Roadnew(t,:)=newp2;
end     
end
Road1=Roadnew;
%%%%%%%%%计算交叉操作之后种群的适应度值%%%%%%%
Dist=zeros(SizeScale,1);
for i=1:SizeScale
     for j=1:(n-1)
         Dist(i)=Dist(i)+DistMatrix(Road1(i,j),Road1(i,j+1));
     end    
     Dist(i)=Dist(i)+DistMatrix(Road1(i,1),Road1(i,n));%计算里程，更新里程矩阵
end 
[MinRoad2 B]=min(Dist);
MaxRoad=max(Dist(:,1));
for i=1:SizeScale
    fitmatrix(i)=fitness(MinRoad2,MaxRoad,Dist(i));%计算适应度函数
end
%%%%%%变异操作（单点变异）%%%%%%%
k=1;
[d p]=sort(fitmatrix(:,1));
while k<=SizeScale
c=randi([1 n],2,1);
pos1(1,:)=c(1,:);%随机产生变异点
pos2(1,:)=c(2,:);%随机产生变异点
rk=rand();
if rk<=pm&&k~=p(SizeScale)%判断是否进行变异
temp=Road1(p(k),pos1);%进行变异操作，更新种群
Road1(p(k),pos1)=Road1(p(k),pos2); 
Road1(p(k),pos2)=temp;
end
k=k+1;
end
%%%%%%%%%计算变异操作后种群的适应度值%%%%%%%%%%%%%
Dist=zeros(SizeScale,1);
for i=1:SizeScale
    for j=1:(n-1)
        Dist(i)=Dist(i)+DistMatrix(Road1(i,j),Road1(i,j+1));
    end
     Dist(i)=Dist(i)+DistMatrix(Road1(i,1),Road1(i,n));
end 
[MinRoad3,C]=min(Dist);
%%%%%%%%%%搜索每一代中的最优路径%%%%%%%%%
if (MinRoad2>MinRoad)&&(MinRoad3>MinRoad)
   MinestRoad=MinRoad;
   D=A;
   Road(D,:)=Road(A,:);
else
    if (MinRoad>MinRoad2)&&(MinRoad3>MinRoad2)
       MinestRoad=MinRoad2;
       D=B;
       Road(D,:)=Roadnew(B,:);
    else
        MinestRoad=MinRoad3;
        D=C;
        Road(D,:)=Road1(C,:);
    end
end
    MinestRoad_fval(iter,1)=MinestRoad;%本代最小里程值
    MinestRoad_opt(iter,:)=Road(D,:);%本代最优路径
    iter=iter+1;
    Road=Road1;
end
[MinestRoad a]=min(MinestRoad_fval);%取路径里程最小值
opt=MinestRoad_opt(a,:)%输出最优路径
fval=MinestRoad%输出最短里程
A=a%得到最优路径的迭代次数
 %%%%%%%%%%绘图%%%%%%%%%
plot(Map(:,1),Map(:,2),'*');
hold on
for c = 1:n
text(Map(c, 1), Map(c, 2), [' ' num2str(c)], 'Color', 'k', 'FontWeight', 'b');
end
XX=Map(MinestRoad_opt(a,:),1);
XX=[XX;Map(MinestRoad_opt(a,1),1)];
YY=Map(MinestRoad_opt(a,:),2);
YY=[YY;Map(MinestRoad_opt(a,1),2)];
plot(XX,YY);
legend('城市','最优路径');
end