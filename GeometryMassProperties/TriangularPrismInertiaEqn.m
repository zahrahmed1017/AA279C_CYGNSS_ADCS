syms a b z h m p

p  = (2*m)/(a*h*b);
dm = -((p*a*b)/h) * (z-h);
x  = ((-z*a)/h) + a;

dIx = ((dm * b^2)/12) + (dm * (z^2 + (b/2)^2));
Ix  = int(dIx, z, 0, h);

dIy = ((dm * x^2)/12) + (dm * (z^2 + (x/2)^2));
Iy  = int(dIy, z, 0, h);

dIz = (dm * (x^2 + b^2)/12) + (dm * ((b/2)^2 + (x/2)^2));
Iz  = int(dIz, z, 0, h);
 
pretty(simplify(Ix));
pretty(simplify(Iy));
pretty(simplify(Iz));
