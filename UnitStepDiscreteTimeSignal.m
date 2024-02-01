clc;
clear all;
close all;

N= input('Enter the range: ');
n= -N:1:N;
y= [zeros(1,N),1,ones(1,N)];
stem(n,y);
grid on;
axis([-(N+1) N+1 -0.5 1.5]); % [-x x -y y]
xlabel('Time');
ylabel('Amplitude of Y');
title('Generating Unit Step Function'); 