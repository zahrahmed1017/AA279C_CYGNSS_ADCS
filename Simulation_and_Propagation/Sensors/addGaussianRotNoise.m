function vec = addGaussianRotNoise(vec, sigma)
    rot_ax = rand(3,1);
    rot_ax = rot_ax ./ norm(rot_ax);

    R_noise = axisAngle2dcm(rot_ax, randn(1)*sigma);

    vec = R_noise * vec;


end