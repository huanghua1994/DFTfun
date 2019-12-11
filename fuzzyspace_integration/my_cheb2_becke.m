function [rad_r, rad_w] = my_cheb2_becke(nCheby, P)
    rad_r = zeros(nCheby, 1);
    rad_w = zeros(nCheby, 1);
    m = nCheby - 1;
    for i = 0 : m
        radx = cos(i * pi / m);
        radr = P * (1 + radx) / (1 - radx);
        radw = (2 * pi / m) * (P^3) * (1 + radx).^2.5 / ((1 - radx).^3.5);
        rad_r(i+1) = radr;
        rad_w(i+1) = radw;
    end
    rad_r = rad_r(end : -1 : 2);
    rad_w = rad_w(end : -1 : 2);
end