function cstr = newlcs(sa,sb,mode)

isseq=false;  % substring
cstr=[];
if nargin<2
    error('MATLAB:lcs','Not enough input arguments.');
elseif nargin==3
    isseq=mode;
end

if isempty(sa) || isempty(sb), return; end

lena=length(sa);
lenb=length(sb);

c = zeros(size(sb));
if ~isseq
    endp=0;
    len=0;
    for i=1:lena
        
        for j=lenb:-1:1
            
            if sa(i) == sb(j)
                
                if(i==1||j==1)
                    c(j)=1;
                else
                    
                    c(j)=c(j-1)+1;
                end
            else
                c(j) = 0;
            end
            
            if c(j) > len
                
                len=c(j);
                endp=j;
            end
        end
    end
    start=endp-len+1;
    cstr=sb(start:endp);
else
    r=zeros(lena+1,lenb+1);
    d=zeros(lena+1,lenb+1);
    for i=2:lena+1
        for j=2:lenb+1
            if sa(i-1) == sb(j-1)
                r(i,j)=1;
                d(i,j)=d(i-1,j-1)+1;
            else
                if d(i-1,j)>=d(i,j-1)
                    d(i,j) = d(i-1,j);
                    r(i,j)=2;
                else
                    d(i,j) = d(i,j-1);
                    r(i,j)=3;
                end
            end
        end
    end
    idx=[];
    while i~= 1 && j~= 1
        switch r(i,j)
            case 1
                idx=[i,idx];
                i=i-1;j=j-1;
                
            case 2
                i=i-1;
            case 3
                j=j-1;
        end
    end
    cstr=sa(idx-1);
    
end