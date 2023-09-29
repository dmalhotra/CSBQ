%Read Matrix
function M=read_mat(fname)
    M=[];
    f1=-1;
    f1=fopen(fname,'r');
    if f1<0
        return;
    end

    n2=2;
    [n1 ,c1]=fread(f1,n2,'uint64');
    if c1~=prod(n2)
        fclose(f1);
        return;
    end

    n1=n1(end:-1:1);
    [M ,c3]=fread(f1,n1','double');
    fclose(f1);
    if c3~=prod(n1)
        M=[];
    end
end
