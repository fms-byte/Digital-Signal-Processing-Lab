clc;
close all;

%orginal signal
n=-6:6;
y=[0 0 0 0 -1 -2 0 2 1 0 0 0 0];
subplot(3,1,1);
stem(n,y);
grid on;
axis([-8 8 -2 2]);
xlabel('n');
ylabel('amplitude ');
title('x[n]');

%shifting the given signal
n1=n-3;
subplot(3,1,2);
stem(n1,y);
grid on;
axis([-8 8 -2 2]);
xlabel('n');
ylabel('amplitude');
title('x[n+3]');

index=1;
for i=1:length(n1)
    if(rem(n1(i), 2)== 0)
        x1(index)=n1(i) / 2;
        y1(index)=y(i);
        index=index+1;
    end

end

%final signal
subplot(3,1,3);
stem(x1,y1);
grid on;
axis([-8 8 -2 2]);
xlabel('n');
ylabel('amplitude');
title('y[n] = x[2n+3]');