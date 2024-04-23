function vx = crossMatrix(v)
% v is a 3x1 column matrix
% outputs the [vx] matrix

vx = [    0,  -v(3),  v(2);
       v(3),      0, -v(1);
      -v(2),   v(1),     0];

end