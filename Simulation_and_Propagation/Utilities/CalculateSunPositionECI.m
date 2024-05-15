function [earth2sun_unitvec_eci, earth2sun_eci] = CalculateSunPositionECI(caldate)

% % This is a simplified (i.e. jank) earth 2 sun position calculator

% HCI to ECI by accounting for obliquity
R_hci2eci = [1,     0,            0;...
             0, cosd(23.45),  sind(23.45);...
             0, -sind(23.45), cosd(23.45)];

sun2earth_hci = OE2HCI(3, juliandate(caldate));
earth2sun_hci = -sun2earth_hci(1:3);
earth2sun_eci = R_hci2eci * earth2sun_hci;

earth2sun_unitvec_eci = earth2sun_eci / norm(earth2sun_eci);

end