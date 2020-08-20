phi = -1.2;
y0 = 1;

dx = 0.1;
x = 0:dx:10;

y = y0*exp(phi*x);

figure;
plot(x,y);

filename = 'einfach.txt';
fid = fopen(filename,'w');
if (fid == -1)
    disp(['unable to open file ' filename]);
    return;
end;

for i = 1:length(x)
    fprintf(fid,'%f %f %f \n',x(i),y(i),1);
end;
fclose(fid);