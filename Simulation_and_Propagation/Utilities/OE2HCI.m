function rv_HCI = OE2HCI(planetID, JD)

    % Step 1: Obtain planetary orbital elements and rate changes:
    [elements, rates] = AA279j2000_planetary_elements(planetID);

    % Step 2: Calculate orbital elements at epoch:
    % oe_epoch = [ a [AU], e, i [deg], OM [deg], om_bar [deg], L [deg] ]
    oe_epoch = zeros(1,6);
    for i = 1:length(elements)
        oe_epoch(i) = elements(i) + (rates(i) * ((JD - 2451545)/36525));
    end

    % Step 3: Convert om_bar to w and L to M:
    OM          = oe_epoch(4);
    om_bar      = oe_epoch(5);
    L           = oe_epoch(6);
    oe_epoch(5) = om_bar - OM; % w
    oe_epoch(6) = L - om_bar;  % M
    

    % Step 4: Convert units:
    oe_epoch(1) = oe_epoch(1) * 149597870.7; % AU to km
    for j = 3:length(oe_epoch) % deg to rad
        oe_epoch(j) = deg2rad(oe_epoch(j));
    end

    % Step 5: Convert Mean anomaly to true anomaly:
    E = M2E(oe_epoch(6), oe_epoch(2), 1e-8);
    true_anom = E2T(E,oe_epoch(2));
    oe_epoch(6) = true_anom;

    % Step 6: Convert oe to HCI position and velocity. Note: Using the
    % OE2ECI function but because our orbital elements are defined wrt the
    % ecliptic and solar system barycenter, the resulting position and vel 
    % will be in HCI not ECI
    rv_HCI = OE2ECI(oe_epoch', 1.3271244004193938e11);

end