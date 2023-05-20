function dxdt=dyn_s_o(t,x,constants)

a = constants(1);
b = constants(2);
c = constants(3);
d = constants(4);

dxdt(1) = a*x(1)*(1-x(1))-b*x(1)*x(2);
dxdt(2) = c*x(2)*(1-x(2))-d*x(1)*x(2);

dxdt = dxdt';

end