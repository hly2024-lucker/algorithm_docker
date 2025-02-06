function Chrom = InitPop(NIND,N,Distance)
%*初始化种群
%输入：
%NIND：种群大小
% N：个体染色体长度(这里为城市的个数)
%输出：
%初始种群
% Chrom = zeros(NIND,N); 	% 用于存储种群
% allow_index = ~ismember(citys_index,tabu);
% allow = citys_index(allow_index);   %待访问的城市集合

for i = 1:N
    City(1,i)=i;    %城市坐标
    %D(i,i) = inf;
end
D = Distance;
Init = zeros(N,4);  %基因库矩阵 
       

for i = 1:N
        [D(i,:),id]=sort(D(i,:));        %city_id:排序后的与i城相邻的城市编号
        for j = 1:4   %找到最邻近的城市距离和坐标
       
            Init(i,j) = id(j);
        end
end         %基因片段建立完成
        
            
    for i = 1:NIND   %第i个个体
        
        index = zeros(1,N); %禁忌表（访问过的城市集合）
       
        total = 1;  %访问过的城市数目 
        s=randi([1,N],1,1);  %随机确定起点城市(1~N)
        index(total) =  s;   %起点城市加入禁忌表 
        Chrom(i,1) = s;
        
        fprintf('begin : %d 已加入禁忌表\n',s);
        allow = City(~ismember(City,index));        %允许访问城市集合
        allow(find(allow==0)) = []; 
        fprintf('allow is \n');
        disp(allow);
        fprintf('禁忌表 is \n');
        disp(index);
        while total<N
            if total == 1
                t = Init(s,2);
                total = total+1;
                index(total) = t; %已访问
                fprintf('s is %d , t is %d , 1st: %d 已加入禁忌表\n',s,t,t);

                Chrom(i,total) = t;
                allow = City(~ismember(City,index));        %更新允许访问城市集合
                allow(find(allow==0)) = []; 
                fprintf('allow is \n');
                disp(allow);
            else
                a = index(total);
                fprintf('total is %d , a is %d\n',total,a);
                
                if ismember(Init(a,1:4),allow)
                    
                    tmp = Init(a,ismember(a,1:4),allow);
                    s2 = tmp(1);        %城市index(total)：最后走过的城市。s2是与它最近的城市
                
                    if ismember(Init(index(1),:),allow)     %左右同时可扩展
                        tmp = Init(index(total),ismember(Init(index(total),:),allow));
                        s1 = tmp(1);        %城市index(1)：城市链最左边。s1是与它最近的城市
                        if s1<=s2    %扩展左边
                            t = s1;
                            index(1,2:total+1) = index(1,1:total);
                            Chrom(i,2:total+1) = Chrom(i,1:total);
                            total = total+1;
                            index(1) = t;
                            fprintf('%d 已加入禁忌表\n',t);
                            Chrom(i,1) = t;
                            allow = City(~ismember(City,index));        %更新允许访问城市集合
                            allow(find(allow==0)) = []; 
                            pause(8);
                        else
                            t = s2;     %仅尾部可拓展
                            total = total+1;
                            index(total) = t; %已访问
                            Chrom(i,total) = t;
                            allow = City(~ismember(City,index));        %更新允许访问城市集合
                            allow(find(allow==0)) = []; 
                        end
                    
                    else
                        t = s2;     %仅尾部可拓展
                        total = total+1;
                        index(total) = t; %已访问
                        Chrom(i,total) = t;
                        allow = City(~ismember(City,index));        %更新允许访问城市集合
                        allow(find(allow==0)) = []; 
                    end
                elseif ismember(Init(index(1),:),allow) %仅头部可扩展
                        t = s1;
                        index(1,2:total+1) = index(1,1:total);
                        Chrom(i,2:total+1) = Chrom(i,1:total);
                        total = total+1;
                        index(1) = t;
                        fprintf('%d 已加入禁忌表\n',t);
                        pause(8);
                        allow = City(~ismember(City,index));        %更新允许访问城市集合
                        allow(find(allow==0)) = []; 
                else
                    [m,n] = size(allow);
                    a = ceil(n*rand(1,1));
                    t = allow(ceil(a));
                    total = total+1;
                    index(total) = t; %已访问
                    %fprintf('%d 已加入禁忌表\n',t);
                    Chrom(i,total) = t;
                     
                    allow = City(~ismember(City,index));        %更新允许访问城市集合
                    allow(find(allow==0)) = []; 
                end
                
            end
                fprintf('allow is \n');
                disp(allow); 
                fprintf('禁忌表 is \n');
                disp(index);
                            
            s = t;
            allow = City(~ismember(City,index));        %更新允许访问城市集合
            allow(find(allow==0)) = []; 
        end  
    end
    for i = 1:NIND    
        p = OutputPath(Chrom(i,:));
    end
    length = zeros(1,NIND);
    for i=1:NIND
        length(i) = PathLength(Distance,Chrom(i,:));
    end
    [m,n] = sort(length);
    tmp = Chrom;
    for i = 1:NIND
        Chrom(n(i),:) = tmp(i,:);
    end
    
    
end


