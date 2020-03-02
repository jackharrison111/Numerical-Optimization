clear all
T = [1,2,3,4];
theta = [2.3,5.1,7.2,9.5];
plot(T,theta);
J=zeros(4);
a = linspace(2.2,2.6,4);
for i = 1:length(T)
  J[i] = (a*T - theta).^2;
endfor

plot(a,J);