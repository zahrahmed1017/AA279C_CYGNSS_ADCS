syms a b z h m p

p  = (2*m)/(a*h*b); % density
dm = -((p*a*b)/h) * (z-h);
x  = ((-z*a)/h) + (2/3)*a;

% set CM at origin and do integrals w/ limits reflecting that

dIx = ((dm * b^2)/12) + (dm * (z^2 )); % use distance from slice CM to prism CM
Ix  = int(dIx, z, -h/3, 2*h/3);

dIy = ((dm * x^2)/12) + (dm * (z^2 + (x/2)^2));
Iy  = int(dIy, z, -h/3, 2*h/3);

dIz = (dm * (x^2 + b^2)/12) + (dm * (  (x/2)^2));
Iz  = int(dIz, z, -h/3, 2*h/3);
 
pretty(simplify(Ix));
pretty(simplify(Iy));
pretty(simplify(Iz));
