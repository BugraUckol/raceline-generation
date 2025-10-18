function [x1] = RK4(x1,omega,f1,h)
  
  k1x1 = f1(x1, omega);
  
  k2x1 = f1(x1+0.5*k1x1*h, omega);
  
  k3x1 = f1(x1+0.5*k2x1*h, omega);
  
  k4x1 = f1(x1+k3x1*h, omega);
   
  x1 = x1 + 1/6*(k1x1+2*k2x1+2*k3x1+k4x1)*h;

end
