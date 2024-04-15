function T = RotationalKineticEnergy(I, w)

T = (1/2)*(I(1,1)*w(1)^2 + I(2,2)*w(2)^2 + I(3,3)*w(3)^2 );

end