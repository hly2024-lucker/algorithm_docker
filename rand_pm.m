function pm = rand_pm(SelCh,Distance,per)	% Distance矩阵的路径长度
NIND = size(SelCh,1);
D = zeros(NIND,NIND);
A = 2;
pm1 = 0.05;
pm2 = 0.002;
%% 计算选择个体的适应度值

ObjV = PathLength(Distance,SelCh);     % 计算路径长度

FitnV = Fitness(ObjV);
f_max = max(FitnV);
f_ave = seek_ave(SelCh,FitnV);
[TobjV,index]=sort(ObjV);

for i = 1:NIND-1
    for j = i+1:NIND
        for n = 1:length(SelCh(i,:))
        	D(i,j) = D(i,j)+abs(SelCh(i,n)-SelCh(j,n));
        end
        D(j,i) = D(i,j);
    end
end
j = index(1);
a = max(max(D));       %% 最大距离的计算

H = mean(D(:));

f = FitnV(per);
    if f>=f_ave
        pm = pm1-(pm1-pm2)/(f_max-f_ave)*A*(f_max-f)*abs(D(i,j)-H)/(a-H);
    elseif f<f_ave
            pm = pm1*A*(D(i,j)-H)/(a-H);
    end
   
end

