function CDWnum = CDWfun(myWing,Cly,M,t_c)
    %Eqn. 15, Eqn. 17, Eqn. 18 of Ref15
    function CdWnum = CdW(M,Lambda50, t_c, Cl)
        Mcrfun = @(lambda, t_c, Cl) (0.95 - t_c./cos(lambda) - Cl/(10*cos(lambda)^2)) /cos(lambda) - (0.1/80)^(1/3);
        Mcr = Mcrfun(Lambda50,t_c, Cl);
        if M <= Mcr
            CdWnum = 0;
        else
            CdWnum = 20*(M - Mcr)^4;
        end
    end
    for i = 1:length(myWing.stripy)
        %Local drag coefficient on each strip
        y = myWing.stripy(i);
        CdWArray = CdW(M,myWing.Lambdax_c(0.5,y),t_c,interp1(Cly(1,:),Cly(2,:),y));
        %disp(interp1(Cly(1,:),Cly(2,:),y))
    end
    CDWnum = sum(CdWArray.*myWing.Sc./myWing.SREF);
end