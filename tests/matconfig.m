p.ti0 = 1.9000e-04;
p.ni0 = 1.5180;
p.ni = 1.5180;
p.tg0 = 1.7000e-04;
p.tg = 1.7000e-04;
p.ng0 = 1.5150;
p.ng = 1.5150;
p.ns = 1.330;
p.lambda = 5.5000e-07;
p.M = 1;
p.NA = 1.4500;
p.pixelSize = 0.05e-06;

xp = [0, 0, 0];
nx = 61;

range = 3;
step = 0.05;

z = (-range/2:step:range/2)*1e-6;
[h, dxp, dyp, dzp] = vectorialPSF(xp, z, nx, p);
imagine(h)