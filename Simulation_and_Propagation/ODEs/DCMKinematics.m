function dcm_d = DCMKinematics(dcm, w)

wx = crossMatrix(w);

dcm_d = -wx * dcm;

end