clc;
clear all;
close all;
m1 = input('Enter the value of x-axis in negative side:');
m2 = input('Enter the value of x-axis in positive side:');
n = -m1:1:m2;
d = [zeros(1,m1), 1 ,zeros(1,m2)];
stem(n,d);
axis([-m1, m2, -0.5, 1.5])
xlabel('n');
ylabel('Amplitude of Y');
title('Unit impulse signal');