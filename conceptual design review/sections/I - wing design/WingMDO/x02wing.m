function wing = x02wing(x)
    wing = WingGeometry();
    wing.cr = x(1);
    wing.ck = x(2);
    wing.ct = x(3);
    wing.s = x(4);
    wing.Lambdain50 = x(5);
    wing.Lambdaout50 = x(6);
    wing.yk = x(7);
    wing.twist_max = x(8);
    wing = wing.calcSref();
end