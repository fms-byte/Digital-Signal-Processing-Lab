close all;
clear all;
clc;

start_value = input('Enter the start value: ');%-6
end_value = input('Enter the end value: ');%6

n1 = start_value:end_value;

y=input("Enter the values of signal = "); %[1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1]

index=1;

n2=(2*start_value):(2*end_value);

value = 2;
for i=1:length(n2)
    x1(i)=n2(i);
    if(rem(n2(i),value)==0)
        y1(i)=y(index);
        index=index+1;
    else
        y1(i)=0;
    end
end

subplot(2,1,1);
stem(n1,y,'r');
xlabel("Time");
ylabel("Amplitude");
grid on;
grid minor;
axis([(start_value-1) (end_value+1) -2 2]);
title("Original signal Y[n]=X[n]");


subplot(2,1,2);
stem(x1,y1,'b');
xlabel("Time");
ylabel("Amplitude");
grid on;
grid minor;
axis([(2*start_value-1) (2*end_value+1) -2 2]);
title("Compresion signal Y[n]=X[n/2]");