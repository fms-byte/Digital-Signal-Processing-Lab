clear;
clc;
n = input ('Insert the value of odd n:');
t = 0:.001:1;
sum = 0;
for f = 1:2:n
    w = sin (2 * pi * f * t);
    sum = sum + w;
end
subplot(1,1,1)
plot(t,sum)
grid on;