function SelCh = Recombin(SelCh,Pc)
%%交叉
%Selch  被选择个体
%Pc    交叉概率
%SelCh  交叉后个体

NSel = size(SelCh,1);
for i = 1:2:NSel - mod(NSel,2)
    if Pc>=rand   %交叉概率
        [SelCh(i,:),SelCh(i+1,:)] = intercross(SelCh(i,:),SelCh(i+1,:));
    end
end