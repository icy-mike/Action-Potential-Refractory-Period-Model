function[output] = RK(f,y0,h)
    k1 = f(y0);
    k2 = f(y0+(1/2)*k1*h);
    k3 = f(y0+(1/2)*k2*h);
    k4 = f(y0+k3*h);
    phi = (1/6)*(k1+2*k2+2*k3+k4)*h;
    output = y0 + phi;
end
