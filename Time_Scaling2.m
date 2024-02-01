close all;
clear all;
clc;
n=-6:6;
y=[1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1];
   
index=1;
value1=2;
for i=1:length(n)
    if(rem(n(i),value1)==0)
        x1(index)=n(i)./value1;
        y1(index)=y(i);
        index=index+1;
    end
end
subplot(3,1,1);
stem(n,y);
xlabel("Time domain");
ylabel("Amplitude");
grid on;
axis([-8 8 -0.5 1.5]);
title("Original signal X[n]");

subplot(3,1,2);
stem(x1,y1,'r');
xlabel("Time domain");
ylabel("Amplitude");
grid on;
grid minor;
axis([-8 8 -0.5 1.5]);
title("Compresion signal Y[n] = X[2n]");

index=1;
n2=-12:12;
value2=2;
for i=1:length(n2)
    x1(i)=n2(i);
    if(rem(n2(i),value2)==0)
        y1(i)=y(index);
        index=index+1;
    else
        y1(i)=0;
    end
end
subplot(3,1,3);
stem(x1,y1,'b');
xlabel("Time domain");
ylabel("Amplitude");
grid on;
grid minor;
axis([-13 13 -0.5 1.5]);
title("Expanding signal Y[n] = X[n/2]");