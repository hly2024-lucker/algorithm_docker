function SelCh = original_Mutate(SelCh,Pm)
%基本位变异
[NSel,L] = size(SelCh);
for i = 1:NSel
    if Pm>=rand
        R = randperm(L);
        SelCh(i,R(1:2)) = SelCh(i,R(2:-1:1));
    end
end

