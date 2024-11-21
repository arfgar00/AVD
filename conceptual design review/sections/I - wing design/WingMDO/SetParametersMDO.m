% Airfoil
polarPath =  '../Airfoils/Airfoil polar data';
XYPath =  '../Airfoils/Airfoil XY data';
polarName = "xf-sc20714-il-1000000.csv";
shapeName = "sc20714.dat.txt";
airfoil = Airfoil();
airfoil = airfoil.readPolar(string(fullfile(polarPath, polarName)));
airfoil = airfoil.readShape(string(fullfile(XYPath, shapeName)));
airfoil = airfoil.interpShape(9);

% Cruise
cruise1 = AirCondition();
cruise1.M = 0.83;
cruise1.h = convlength(28000, 'ft','m');

cruise2 = AirCondition();
cruise2.V = convvel(400, 'kts', 'm/s');
cruise2.h = convlength(18000, 'ft','m');
cruise2 = cruise2.calcM();

% Design CL
CLdesC1 = 0.4623;
CLdesC2 = 0.3008;

% Misc
C_cruise = 0.5/3600;
bodyDiameter = 6.38;
Wto0 = 4.2188*1e5*9.81;
Wpayload = 8.7176e3*9.81;
Swet_SrefWing = 2.2;
bodyLength = 72;
%Swet_SrefBody = pi*bodyDiameter*bodyLength*1.1/Srefdes;
l_d = bodyLength/bodyDiameter;
Srefdes = 548.7;
W_Sdes = Wto0/Srefdes;
Wto = Wto0;
nult = 3.75;

% Weight of Zero Fuel without Wing
Ww0 = 35220.312*9.81; % Guessed from wing0
W_ZF_Nowing = Wto0 - Wto0 * 0.42 - Wpayload - Ww0;

% Wing span
bdes = 64.99;

% For LLE
Damping = 0.04;
tolerance = 1e-2;
