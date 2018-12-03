a=0.1:0.2:0.9;
b=0.1:0.1:0.9;
inum=1;
for i=1: length(a)
    for j=1: length(b)
        saveDir=['G:\TEST\' 'D_' num2str(a(i)) '_T_'  num2str(b(j)) ] ;
        system(['mkdir ' saveDir]);
        cd(saveDir);
        data=[a(i); b(j)];
        save('variable.txt', 'data','-ascii');
    end
end