%y[n]=x[n-3] and z[n]=x[n+2]
clc;
close all;

%orginal signal
n=-3:3;
y=[-2 0 1 -3 2 -1 3];
subplot(3,1,1);
stem(n,y);
grid on;
axis([-6 7 -4 4]);
xlabel('n');
ylabel('amplitude ');
title('x[n]');

%shifting the given signal
n1=n+3;
subplot(3,1,2);
stem(n1,y);
grid on;
axis([-6 7 -4 4]);
xlabel('n');
ylabel('amplitude');
title('y[n] = x[n-3]');

%shifting the given signal
n2=n-2;
subplot(3,1,3);
stem(n2,y);
grid on;
axis([-6 7 -4 4]);
xlabel('n');
ylabel('amplitude');
title('z[n] = x[n+2]');