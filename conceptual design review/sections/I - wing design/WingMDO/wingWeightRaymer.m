function Ww = wingWeightRaymer(mywing,myairfoil,Wto,Nult,Scswin,Scswout)
    %Scsw is the surface area of the control surfaces mounted on the wing,
    %on inner section and outer section respectively
    Wwing = @(Wdg,Nz,Sw,A,lambda,Scsw,Lambda,t_croot) 0.0051*((Wdg*Nz)^0.557*Sw^0.649*A^0.5*(1+lambda)^0.1*Scsw^0.1)/(t_croot^0.4*cos(Lambda));
    ARin = (2*mywing.yk)^2 / mywing.Sin;
    ARout =  (2*(mywing.s-mywing.yk))^2 / mywing.Sout;
    lambdain = mywing.ck/mywing.cr;
    lambdaout = mywing.ct/mywing.ck;
    Wwin = Wwing(Wto,Nult,mywing.Sin,mywing,ARin,lambdain,Scswin,mywing.Lambdain50,myairfoil.t_c)
    Wwout = Wwing(Wto,Nult,mywing.Sout,mywing,ARout,lambdaout,Scswout,mywing.Lambdaout50,myairfoil.t_c)
    Ww = Wwin + Wwout;
end