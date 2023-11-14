function [y]=f(x)
 %y=tanh(x);
% y=x;
y=0.5*(abs(x+1)-abs(x-1));
end