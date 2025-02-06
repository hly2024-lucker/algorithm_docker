function avefit = seek_ave(Chrom,FitnV)

NIND = size(Chrom,1);
cumfit = 0;
if NIND>0
    for i = 1:NIND
        cumfit = cumfit + FitnV(i);
    end
end
avefit = cumfit/NIND;   %平均适应度


end