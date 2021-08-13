count = 1;
file_name = "info2.txt";
data = readtable(file_name);

figure(count)
p = islocalmin(data.cost);
plot(data.psfSize,p);
title("psf vs cost");
count = count +1;
saveas(gcf,"psf vs cost.png");

figure(count)
p1 = islocalmin(data.outVarlap1);
p2 = islocalmin(data.outVarlapnorm1);
plot(data.psfSize,p1,data.psfSize,p2);
title("psf vs lap1")
count = count +1;
saveas(gcf,"psf vs lap1.png");

figure(count)
p1 = islocalmin(data.outMean);
p2 = islocalmin(data.outMeannorm);
plot(data.psfSize,p1,data.psfSize,p2);
title("psf vs mean");
count = count +1;
saveas(gcf,"psf vs mean.png");

figure(count)
p1 = islocalmax(data.outVarlap2);
p2 = islocalmax(data.outVarlapnorm2);
plot(data.psfSize,p1,data.psfSize,p2);
title("psf vs lap2");
count = count +1;
saveas(gcf,"psf vs lap2.png");

disp("----");
disp(file_name);
n = islocalmin(data.cost) & ( surroundone(islocalmin(data.outVarlap1)) | surroundone(islocalmin(data.outMean)) | surroundone(islocalmax(data.outVarlap2)));
[l,p] = islocalmin(data.cost);
[~,idx] = max(p);
n(idx) = 1;
disp("prediction size : ");
disp(data.psfSize(n)');

function b = surroundone(a)
    b = a;
    b(1) = a(1) | a(2);
    for i = 2:numel(a)-1
    b(i) = a(i) | a(i+1) | a(i-1);
    end
end

