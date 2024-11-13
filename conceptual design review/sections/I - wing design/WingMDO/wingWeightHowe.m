%Weight estimation given in Ref19 Eqn. 5-7, Howe
%nult: ultimate load factor
%Wpay: payload weight
%Wto: take off weight
%t_cr: t/c at root
%Wpay, Wto are in N
%VD: maximum structural design velocity
function Ww = wingWeightHowe(mywing, Wpay, Wto, nult, t_cr, VD, range)
    %C1 = 0.00072 - 0.0005*(2.08  + 0.00038*range)*Wpay*1e-6;
    PAX = 516;
    C1 = 0.00072 - 0.0005 * (270  + 0.05*range) * PAX * 1e-6;
    C5 = 1.16;
    %H is the sec Lambda (1+2lambda)/(3+3lambda) in the original formula.
    %This is weighted averaged of two sections.
    lambdain = mywing.ck/mywing.cr;
    lambdaout = mywing.ct/mywing.ck;
    lambda = mywing.ct / mywing.cr;
    g = 9.81;
    %H = mywing.Sin/mywing.SREF*sec(mywing.Lambdain50)*(1+2*lambdain)/(3+ 3*lambdain) + mywing.Sout/mywing.SREF*sec(mywing.Lambdaout50)*(1+2*lambdaout)/(3+ 3*lambdaout);
    %Ww = 9.81*C1/C5*(mywing.AR^0.5*mywing.SREF^1.5*H *Wto/(9.81*mywing.SREF)*nult*(VD/t_cr)^0.5)^0.9;
    Ww = g*C1/C5*(mywing.AR^0.5 * mywing.SREF^1.5 / cos(mywing.Lambdaout50) * (1 + 2*lambda)/(3 + 3*lambda) * Wto/9.81/mywing.SREF * nult * (VD/t_cr)^0.5 )^0.9;
end