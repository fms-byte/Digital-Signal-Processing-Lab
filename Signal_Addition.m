clear all;
clc;
x1=[-5 -4 -3 -2 -1 0];
y1=[2 5 4 6 3 5];
x2=[-2 -1 0 1 2];
y2=[8 9 2 5 6];
 
% Draw the second signal.
subplot(3,1,1);
stem(x1,y1);
grid on;
grid minor;
axis([-10 10 -8 8]);
 
% Draw the second signal.
subplot(3,1,2);
stem(x2,y2);
grid on;
grid minor;
axis([-10 10 -2 16]);
n=min(min(x1),min(x2)):1:max(max(x1),max(x2));
 
% This function is for the addition the two signal .
[y] = add_function(n,x1,x2,y1,y2);
 
% This is for the plot the added signal.
subplot(3,1,3);
stem(n,y);
grid on;
grid minor;
axis([-10 10 -2 16]);


function[y] = add_function(n,x1,x2,y1,y2)
 
m1=zeros(1,length(n));
m2=zeros(1,length(n));
temp=1;
for i=1:length(n)
    if(n(i)>=min(x1) & n(i)<=max(x1))
        m1(i)=y1(temp);
        temp=temp+1;
    else
        m1(i)=0;
    end
end
temp=1;
for i=1:length(n)
    if(n(i)>=min(x2) & n(i)<=max(x2))
        m2(i)=y2(temp);
        temp=temp+1;
    else
        m2(i)=0;
    end
end
 
y=m1+m2;
end
 

