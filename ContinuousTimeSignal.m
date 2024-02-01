% Generation of Continuous time signals 
% 2sin(2πτ-π/2)
T = -5:0.05:5;
x=2*sin((2*pi*T) - (pi/2));
plot(T,x)
grid on;
axis ([-6 6 -3 3])
ylabel ('x(t)')
xlabel ('Time(sec)')
title ('Figure 2.1')