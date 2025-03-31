function theta_wrapped = wrapToPi(theta, is_rad)
% is_rad: 1 if theta is already in [rad]
% 0 if theta is in [deg]
% need to covert the rad angle and make it located in -pi to pi
% return: the wrapped angle in rad.

    if ~is_rad
        theta = theta/180*pi;
    end
    theta_wrapped = mod(theta + pi, 2*pi) - pi;
end