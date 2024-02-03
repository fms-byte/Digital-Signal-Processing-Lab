%n = deg2rad(-180:30:180);
n = -2:1:2;
x = [5,5,4,3,1];
%x = cos(n) + sin(n) + (sin(n) .* cos(n));
%x = 1 + n + 3.*(n.*n) + 5.*(n.*n.*n) + 9.*(n.*n.*n.*n);
%x = 1 + n.*(cos(n)) + (n.*n).*(sin(n)) + (n.*n.*n).*(sin(n));
%x = (1 + n.^3) .* (cos(10.*n)).^3;

disp(x);
% Creating mirrored versions for negative indices
x_mirror = fliplr(x); %x_mirror = [1,4,3,5,5]

% even and odd components
xe = (x + x_mirror) / 2; %xe= ( x(n)+x(-n) )/2;
xo = (x - x_mirror) / 2; %xo= ( x(n)-x(-n) )/2;

% Plotting
subplot(4,1,1);
stem(n, x); 
grid on; 
axis([-3.5 3.5 min(x)-1 max(x)+1]);
xlabel('n'); 
ylabel('Amplitude'); 
title('Original Signal');

subplot(4,1,2);
stem(n, x_mirror); 
grid on; 
axis([-3.5 3.5 min(x_mirror)-0.5 max(x_mirror)+0.5]);
xlabel('n'); 
ylabel('Amplitude'); 
title('Reversed Signal');

subplot(4,1,3); 
stem(n, xe, 'b'); 
grid on; 
axis([-3.5 3.5 min(xe)-1 max(xe)+1]);
xlabel('n'); 
ylabel('Amplitude'); 
title('Even Signal');

subplot(4,1,4); 
stem(n, xo, 'b'); 
grid on; 
axis([-3.5 3.5 min(xo)-1 max(xo)+1]);
xlabel('n'); 
ylabel('Amplitude'); 
title('Odd Signal');
