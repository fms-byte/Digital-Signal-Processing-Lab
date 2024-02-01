clc;
clear all;
close all;
 
frequency =input('Enter the frequency for the signal:\n');
fprintf("Enter the  frequency and it must be greater or less than Nyquist frequency\n");
oversampling=input('');
%fprintf("Enter the UnderSampling frequency and it must be Less than %d\n",2*frequency);
%undersampling=input('');

Time_Period = 1/frequency;
 
tmin=-0.05;
tmax=0.05;

time = linspace(tmin,tmax,400);
amplitude = cos(2*pi*frequency*time);
 
subplot(4,1,1);
plot(time,amplitude);
grid on; grid minor;
xlabel("Time");
ylabel("Amplitude");
title("orginal Signal");

%Nyquist rate Sampling Part.
nyquist_frequency = 2*frequency;
time1=tmin:(1/nyquist_frequency):tmax;
amplitude1 = cos(2*pi*frequency*time1);
subplot(4,1,2);
plot(time,amplitude);
hold on;
grid on; grid minor;
plot(time1,amplitude1);
title("NyQuist Sampling");
hold off;

if oversampling>nyquist_frequency
    %OverSampling Part.
    time1=tmin:(1/oversampling):tmax;
    amplitude1 = cos(2*pi*frequency*time1);
    subplot(4,1,3);
    plot(time,amplitude);
    hold on;
    plot(time1,amplitude1);
    grid on; grid minor;
    xlabel("Time");
    ylabel("Amplitude");
    title("OverSampling");
    hold off;
else
    %OverSampling Part.
    time1=tmin:(1/oversampling):tmax;
    amplitude1 = cos(2*pi*frequency*time1);
    subplot(4,1,3);
    plot(time,amplitude);
    hold on;
    plot(time1,amplitude1);
    grid on; grid minor;
    xlabel("Time");
    ylabel("Amplitude");
    title("UnderSampling");
    hold off;
end
