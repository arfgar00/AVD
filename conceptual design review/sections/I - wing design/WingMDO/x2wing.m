function wing = x2wing(x)
    global bdes
    wing = WingGeometry();
    wing.cr = x(1);
    wing.ck = x(2);
    wing.ct = x(3);
    wing.s = bdes/2;
    wing.N = 501;
    wing.Lambdain50 = x(4);
    wing.Lambdaout50 = x(5);
    wing.yk = x(6);
    wing.twist_max = x(7);
    wing = wing.calcSref();
    wing = wing.createStrips();
    wing = wing.calcSc();
end