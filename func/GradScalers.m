function res = GradScalers(azimuth, polar)

resz = cos(polar);
resy = sin(polar)*sin(azimuth);
resx = sin(polar)*cos(azimuth);

res = [resx, resy, resz];

end