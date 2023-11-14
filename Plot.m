clc;
close all;
clear;
[a,aa]=staticETM;
[b,bb]=dynamicETM;
[c,cc]=periodicETM;
[d,dd]=DETM;

h=0.001;
figure(1);
subplot(2,2,1);
stem(a*h,aa*h,'o','r');
% axis([0 10 0 1.1]);
xlabel('(a)'); 
ylabel('The release time interval');
grid on;
subplot(2,2,2);
stem(b*h,bb*h,'o','g');
% axis([0 10 0 1.8]);
xlabel('(b)'); 
ylabel('The release time interval');
grid on;
subplot(2,2,3);
stem(c*h,cc*h,'o','b');
% axis([0 10 0 0.7]);
xlabel('(c)'); 
ylabel('The release time interval');
grid on;
subplot(2,2,4);
stem(d*h,dd*h,'o','m');
% axis([0 10 0 2]);
xlabel('(d)'); 
ylabel('The release time interval');
grid on;







