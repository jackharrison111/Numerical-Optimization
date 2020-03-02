clear all

function f = Q(x)
  A = [4 2; 2 3];
  b = [3 1];
  c = -1;
  f = 1/2*x'*A*x + b*x + c;
endfunction

function g = grad(x)
  A = [4 2; 2 3];
  b = [3 1];
  g = 1/2*A'*x + 1/2*A*x + b';
endfunction

A = [4 2; 2 3];
b = [3 1];
c = -1;
x_0 = [-1,-1]';
y = Q(x_0);
g_k = grad(x_0);
x_k0 = x_0;
x_k1 = x_k0 - ((g_k'*g_k)/(g_k'*A*g_k))*g_k;


function min = SD(start_guess, tolerance)
  
  x_k0 = start_guess;
  g_k = grad(x_k0);
  A = [4 2; 2 3];
  x_k1 = x_k0 - (g_k'*g_k)/(g_k'*A*g_k)*g_k;
  while abs(x_k1-x_k0) > tolerance
    x_k0 = x_k1;
    g_k = grad(x_k0);
    a_k = (g_k'*g_k)/(g_k'*A*g_k);
    x_k1 = x_k0 - a_k*g_k;
  endwhile
  min = x_k1
endfunction
tic
test = SD(x_k0, 10e-3)
toc

