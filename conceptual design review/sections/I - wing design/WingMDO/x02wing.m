function wing = x02wing(x, cr, s)
    wing = WingGeometry();
    wing.cr = cr;
    wing.ck = x(1);
    wing.ct = x(2);
    wing.s = s;
    wing.Lambdain50 = x(3);
    wing.Lambdaout50 = x(4);
    wing.yk = x(5);
    wing.twist_max = x(6);
    wing = wing.calcSref();
end