dt={'1d-3', '5d-4', '1d-4'};
rtol={'1d-4', '1d-5', '1d-6'}

fid=fopen(sprintf('dt_%s_rtol_%s.dat', dt{end}, rtol{end}), 'r');
a=fscanf(fid, '%f', inf);
ref = reshape(a, 540, length(a)/540);
fclose(fid);
ref2 = [ref ref(:,end)];

nlong = size(ref2, 2);
dr=zeros(nlong/10,length(dt)*length(rtol));
labels={};

n=1;

for i=1:length(dt)
    dn = eval([dt{1} '/' dt{i}]);
    for j=1:length(rtol)
        fid=fopen(sprintf('dt_%s_rtol_%s.dat', dt{i}, rtol{j}), 'r');
        a=fscanf(fid, '%f', inf);
        r = reshape(a, 540, length(a)/540);
        fclose(fid);
        if (size(r, 2) == size(ref, 2))
            dr_ = r-ref;
            dr_ = [dr_ dr_(:,end)];

            dr(:,n) = sqrt(sum(dr_(:,1:10:end).^2, 1)/180);
        else
            dr(:,n) = sqrt(sum((r(:,1:dn:end)-ref2(:,1:10:end)).^2, 1)/180);
        end
        labels{n} = sprintf('dt=%s, rtol=%s', dt{i}, rtol{j});
        n = n+1;
    end
end
        