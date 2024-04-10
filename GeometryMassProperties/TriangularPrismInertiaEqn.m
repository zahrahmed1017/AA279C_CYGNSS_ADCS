syms a b z h m

p  = (2*m)/(a*h*b); % density
x  = ((-z*a)/h) + (a/3);
dm = p * (x + (a/3)) * b;

% set CM at origin and do integrals w/ limits reflecting that

dIx = ((dm * b^2)/12) + (dm * z^2); % use distance from slice CM to prism CM
Ix  = int(dIx, z, -h/3, 2*h/3);

dIy = ((dm * (x + (a/3))^2)/12) + (dm * (z^2 + ((a/6) - (x/2))^2));
Iy  = int(dIy, z, -h/3, 2*h/3);

dIz = (dm * ((x + a/3)^2 + b^2)/12) + (dm * ((a/6) - (x/2))^2);
Iz  = int(dIz, z, -h/3, 2*h/3);
 
pretty(simplify(Ix));
pretty(simplify(Iy));
pretty(simplify(Iz));
