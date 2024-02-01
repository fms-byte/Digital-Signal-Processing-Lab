clc;
clear all;
m1 = input('Enter the value of x-axis in negative side:');
m2 = input('Enter the value of x-axis in positive side:');
n = m1:m2;
x = (n==0);%it works as if statement like n=-5:5( 0 0 0 0 0 1 0 0 0 0 0 0)
stem(n,x);
xlabel('n');
ylabel('amplitude');
title('Unit impulse signal');