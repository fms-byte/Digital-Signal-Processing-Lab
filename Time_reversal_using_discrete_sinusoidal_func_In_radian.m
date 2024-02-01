%close all
%clc
 
t1 = 0:0.2:2*pi; %values of x-axis in radian
x1 = sin(t1);    %values of y-axis
x2 = fliplr(x1); %fliplr() -> this function gives the flipped result; 
                 %lr means left right ...flipud() ud means up down 
t2 = -fliplr(t1); % time values must be flipped and negated
subplot(2,1,1)
stem(t1,x1,'LineWidth',2)
xlim([-10 10])
title('Original Signal') %\bf\fontsize{25}
xlabel('Samples (t)')
ylabel('Amplitude (sin(t))')
grid on;
ax = gca;
ax.XAxis.FontSize = 15;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 15;
ax.YAxis.FontWeight = 'bold';
 
subplot(2,1,2)
stem(t2,x2,'LineWidth',2)
xlim([-10 10])
title('Time Reversed Signal')
xlabel('Samples')
ylabel('Amplitude')
grid on;
ax = gca;
ax.XAxis.FontSize = 15;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 15;
ax.YAxis.FontWeight = 'bold';