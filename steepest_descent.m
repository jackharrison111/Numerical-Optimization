function min = steepest_descent(@func, @grad, start_guess, tolerance, learning_rate):
  x_k0 = start_guess;
  g_k = grad(x_k0);
  a_k = learning_rate;
  x_k1 = x_k0 - learning_rate*g_k;
  while abs(x_k1-x_k0) > tolerance
    x_k0 = x_k1;
    g_k = grad(x_k0)
    x_k1 = x_k0 - a_k*g_k;
  endwhile
  min = x_k1
endfunction
