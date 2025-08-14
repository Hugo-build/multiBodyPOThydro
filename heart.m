
clear all
clc
f = @(x,a)real(abs(x).^(2/3)+exp(1)/3.*(pi-x.^2).^0.5.*sin(a*pi.*x));
x=-2:0.01:2;
for a=0:1:50
    y=f(x,a);
    plot(x,y,'r','LineWidth',2)
    grid on
    axis([-2 2 -1.5 2.5]);
    pause(0.1)
end


