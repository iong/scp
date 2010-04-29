t=1:1040;
r=zeros(1040,18);
for i=1:1040
    c=load(sprintf('%06.0f.dat',i));
    r(i,1:9) = c(29,1:9);
    r(i,10:18) = c(30,1:9);
    %scatter3(c(:,1), c(:,2), c(:,3));
    %title(sprintf('i = %d', i))
    %pause
end