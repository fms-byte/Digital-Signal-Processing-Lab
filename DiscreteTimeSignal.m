% Generation of discrete time signals 
n = -5:5;
x = [0 0 1 1 -1 0 2 -2 3 0 -1];
stem (n,x);
axis ([-6 6 -3 3]);
xlabel ('n'); 
ylabel ('x[n]'); 
title ('Figure 2.2');
  