x = 0:0.1:10;
y = sin (x);
z = cos (x);
subplot (3,1,1);
plot (x,y,'b');
grid on;
xlabel("x");
ylabel("y=sin(x)");
subplot (3,1,2);
plot (x,z,'r');
grid on; 
xlabel("x");
ylabel("z=cos(x)");
subplot (3,1,3);
stem (x,z,'r');
grid on;
hold on;
subplot (3,1,3);
stem (x,y,'b');
xlabel("x");
ylabel("y=sin(x), z=cos(x)");