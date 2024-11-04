%Ref19, Eqn. 5-13, LTH

function Ww = wingWeightLTH(mywing,t_c,Wto)
    H = mywing.Sin/mywing.SREF/cos(mywing.Lambdain50) + mywing.Sout/mywing.SREF/cos(mywing.Lambdaout50);
    Ww = 2.20013*9.81*1e-4*(401.146*mywing.SREF^1.31 + (Wto / 9.81)^1.1038) * (t_c)^(-0.5)*mywing.AR^1.5*H;
end