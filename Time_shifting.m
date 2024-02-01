clc;
clear;
 
n = -5:5;
x = [0 -1 -.5 .5 1 1 1 1 .5 0 0];
subplot(3,1,1);
stem (n,x);
xlabel('Time Sample');
ylabel('Amplitude');
title('Original Signal');
axis([-7 7 min(x)-2 max(x)+2]);
grid on;
grid minor;
 
m = n+4; 
subplot(3,1,2);
stem (m,x);
xlabel('Time Sample');
ylabel('Amplitude');
title('Time right shifted signal');
axis([-7-2+4 7+2+4 min(x)-2 max(x)+2]);
grid on;
grid minor;

l = n-4; 
subplot(3,1,3);
stem (l,x);
xlabel('Time Sample');

ylabel('Amplitude');
title('Time left shifted signal');
axis([-7-2-4 7+2-4 min(x)-2 max(x)+2]);
grid on;
grid minor;