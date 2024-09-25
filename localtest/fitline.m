x = [0.1 0.2 0.5 1 2 3];
y = [4.06 4.04  3.266 2.53 2.098 2.08];

%
% x =minOverlap:minOverlap+length(parameters1)-1;
cnu = polyfit(x,y,1);
y_est1 = polyval(cnu,x);

figure
plot(x,y)
hold on
plot(x,y_est1,'r--','LineWidth',2)