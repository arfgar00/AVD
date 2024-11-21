function wing = x2wing(x)
    bdes = 65;
    wing = WingGeometry();
    wing.cr = x(1);
    wing.ck = x(2);
    wing.ct = x(3);
    wing.s = bdes/2;
    wing.N = 301;
    wing.Lambdaout50 = x(4);
    wing.yk = x(5);
    wing.twist_max = x(6);
    wing.Lambdain50 = wing.returnModedLambda();
    wing = wing.calcSref();
    wing = wing.createStrips();
    wing = wing.calcSc();
end