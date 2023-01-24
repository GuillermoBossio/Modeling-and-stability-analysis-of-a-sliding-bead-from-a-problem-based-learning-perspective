function [SOL] = ODEIntegration (SFrenet,Solver,Parameter)

Options = odeset ('AbsTol', Solver.TOL);

if strcmpi(Solver.IFLAG, 'ForwardEuler') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    Zout(1:2,1) = CI;
    for i = 2:length(Parameter.Time)
        DZ = Fcn (Parameter.Time(i-1), Zout(:,i-1),Parameter, SFrenet);
        Zout(:,i) = Zout(:,i-1) + DZ * Solver.DT;
    end
    Zout = transpose(Zout);
    SOL.Tout = Parameter.Time;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'BackwardEuler')

    AUX = EtaDerivatives (SFrenet);
    SOL.Z(1:2,1) = Parameter.CI;

    TOL = Solver.TOL; Iter = Solver.Iter; DT = Solver.DT; g = Parameter.g;

    for i = 2:length(Parameter.Time)

        Error = 1;
        k = 1;
        Z1 = SOL.Z(1:2,i-1);

        while Error >= TOL && k <= Iter

            b  = SFrenet.Funb (Z1(1,1));
            Db = SFrenet.FunDb (Z1(1,1));
            DDb = AUX.FunDDb (Z1(1,1));
            Dinb = AUX.FunDinvb (Z1(1,1));
            Dx3 = SFrenet.FunDR (Z1(1,1)); Dx3 = Dx3(3,1);
            DDx3 = AUX.FunDDx3 (Z1(1,1));
            Dinb2 = AUX.FunDinvb2 (Z1(1,1));

            F = [Z1(2,1); -Db/b * Z1(2,1)^2 - g*Dx3/b^2];

            RHS = -1 * ( Z1 - SOL.Z(:,i-1) - F*DT);

            JJ      = zeros(2,2);
            JJ(1,1) = 1;
            JJ(1,2) = -DT;
            JJ(2,1) = (DDb/b*Z1(2,1)^2 + Db*Dinb*Z1(2,1)^2 + g*DDx3/b^2 + g*Dx3*Dinb2)*DT;
            JJ(2,2) = 1 + 2*Db/b*Z1(2,1)*DT;

            SOLN = JJ \ RHS;

            Error = norm(SOLN,"inf");

            Z1 = SOLN + Z1;

            k = k + 1;

        end

        SOL.Z (1:2,i)  = Z1;
        SOL.Iter(i)    = k-1;
        SOL.Error(i)   = Error;

    end

    SOL.Z = transpose (SOL.Z);
    SOL.Tout = Parameter.Time;

elseif strcmpi(Solver.IFLAG, 'RungeKutta2') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    Zout(1:2,1) = CI;
    for i = 2:length(Parameter.Time)
        K1 = Fcn (Parameter.Time(i-1), Zout(:,i-1),Parameter, SFrenet);
        Arg2 = Zout(:,i-1) + Solver.DT*K1;
        K2 = Fcn (Parameter.Time(i), Arg2, Parameter, SFrenet);
        Zout(:,i) = Zout(:,i-1) + Solver.DT* (0.5*K1 + 0.5*K2);
    end
    Zout = transpose(Zout);
    SOL.Tout = Parameter.Time;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'RungeKutta3') 

    CI    = Parameter.CI; DT = Solver.DT;
    Fcn   = str2func(Solver.FName);
    Zout(1:2,1) = CI;
    for i = 2:length(Parameter.Time)
        K1 = Fcn (Parameter.Time(i-1), Zout(:,i-1),Parameter, SFrenet);
        Arg2 = Zout(:,i-1) + DT*K1/2;
        TT   = Parameter.Time(i-1) + DT/2;
        K2 = Fcn (TT, Arg2, Parameter, SFrenet);
        Arg2 = Zout(:,i-1) - DT*K1 + 2*DT*K2;
        K3 = Fcn (Parameter.Time(i), Arg2, Parameter, SFrenet);
        Zout(:,i) = Zout(:,i-1) + DT/6 * (K1 + 4*K2 + K3);
    end
    Zout = transpose(Zout);
    SOL.Tout = Parameter.Time;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'RungeKutta4') 

    CI    = Parameter.CI; DT = Solver.DT;
    Fcn   = str2func(Solver.FName);
    Zout(1:2,1) = CI;
    for i = 2:length(Parameter.Time)
        K1 = Fcn (Parameter.Time(i-1), Zout(:,i-1),Parameter, SFrenet);
        Arg2 = Zout(:,i-1) + DT*K1/2;
        TT   = Parameter.Time(i-1) + DT/2;
        K2 = Fcn (TT, Arg2, Parameter, SFrenet);
        Arg2 = Zout(:,i-1) + DT*K2/2;
        K3 = Fcn (TT, Arg2, Parameter, SFrenet);
        Arg2 = Zout(:,i-1) + DT*K3;
        K4 = Fcn (Parameter.Time(i), Arg2, Parameter, SFrenet);
        Zout(:,i) = Zout(:,i-1) + DT * (K1/6 + K2/3 + K3/3 + K4/6);
    end
    Zout = transpose(Zout);
    SOL.Tout = Parameter.Time;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'AdamsBashforth2') 

    CI    = Parameter.CI; DT = Solver.DT;
    Fcn   = str2func(Solver.FName);
    Zout(1:2,1) = CI;
    
    K1 = Fcn (Parameter.Time(1), CI,Parameter, SFrenet);
    Arg2 = Zout(:,1) + Solver.DT*K1;
    K2 = Fcn (Parameter.Time(2), Arg2, Parameter, SFrenet);
    Zout(:,2) = Zout(:,1) + Solver.DT* (0.5*K1 + 0.5*K2);    
    
    for i = 3:length(Parameter.Time)
        F1 = Fcn (Parameter.Time(i-2), Zout(:,i-2),Parameter, SFrenet);
        F2 = Fcn (Parameter.Time(i-1), Zout(:,i-1),Parameter, SFrenet);
        Zout(:,i) = Zout(:,i-1) + DT/2 * (3*F2 - F1);
    end
    Zout = transpose(Zout);
    SOL.Tout = Parameter.Time;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'ODE15sMatlab') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    TSPAN = [0 Solver.SimT];

    [Tout, Zout] = ode15s (Fcn, TSPAN, CI, Options, Parameter, SFrenet);

    SOL.Tout = Tout;
    SOL.Z    = Zout;
     
elseif strcmpi(Solver.IFLAG, 'ODE23Matlab') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    TSPAN = [0 Solver.SimT];

    [Tout, Zout] = ode23 (Fcn, TSPAN, CI, Options, Parameter, SFrenet);

    SOL.Tout = Tout;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'ODE23sMatlab') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    TSPAN = [0 Solver.SimT];

    [Tout, Zout] = ode23s (Fcn, TSPAN, CI, Options, Parameter, SFrenet);

    SOL.Tout = Tout;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'ODE23tMatlab') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    TSPAN = [0 Solver.SimT];

    [Tout, Zout] = ode23t (Fcn, TSPAN, CI, Options, Parameter, SFrenet);

    SOL.Tout = Tout;
    SOL.Z    = Zout;    

elseif strcmpi(Solver.IFLAG, 'ODE45Matlab') 

    CI    = Parameter.CI;
    Fcn   = str2func(Solver.FName);
    TSPAN = [0 Solver.SimT];

    [Tout, Zout] = ode45 (Fcn, TSPAN, CI, Options, Parameter, SFrenet);

    SOL.Tout = Tout;
    SOL.Z    = Zout;

elseif strcmpi(Solver.IFLAG, 'MidPointRule')

    AUX = EtaDerivatives (SFrenet);
    SOL.Z(1:2,1) = Parameter.CI;

    TOL = Solver.TOL; Iter = Solver.Iter; DT = Solver.DT; g = Parameter.g;

    for i = 2:length(Parameter.Time)

        Error = 1;
        k = 1;
        Z1 = SOL.Z(1,i-1);

        Db0 = SFrenet.FunDb (SOL.Z(1,i-1));
        b0  = SFrenet.Funb (SOL.Z(1,i-1));
        Dx30 = SFrenet.FunDR (SOL.Z(1,i-1)); Dx30 = Dx30(3,1);

        while Error >= TOL && k <= Iter

            DZ1 = 2/DT * (Z1 - SOL.Z(1,i-1)) - SOL.Z(2,i-1);
            b  = SFrenet.Funb (Z1);
            Db = SFrenet.FunDb (Z1);
            DDb = AUX.FunDDb (Z1);
            Dinb = AUX.FunDinvb (Z1);
            Dx3 = SFrenet.FunDR (Z1); Dx3 = Dx3(3,1);
            DDx3 = AUX.FunDDx3 (Z1);
            Dinb2 = AUX.FunDinvb2 (Z1);

            RHS = -1 * ( (DZ1 - SOL.Z(2,i-1))/DT + 0.5 * (Db/b*DZ1^2 + Db0/b0*SOL.Z(2,i-1)^2) + ...
                0.5*g * (Dx3/b^2 + Dx30/b0^2) );

            JJ   = (2/DT^2 + (2/DT)*Db/b*DZ1 + 0.5*DDb*DZ1^2/b + 0.5*Dinb*DZ1^2*Db + 0.5*g*DDx3*1/b^2 + 0.5*g*Dx3*Dinb2);

            SOLN = RHS / JJ;

            Error = abs(SOLN);

            Z1 = SOLN + Z1;

            k = k + 1;

        end

        SOL.Z (1,i)  = Z1;
        SOL.Z (2,i)  = 2/DT * (Z1 - SOL.Z(1,i-1)) - SOL.Z(2,i-1);
        SOL.Iter(i)    = k-1;
        SOL.Error(i)   = Error;

    end

    SOL.Z = transpose (SOL.Z);
    SOL.Tout = Parameter.Time;

end

%% Postprocessing 

E3 = [0;0;1];

for i = 1:size(SOL.Z,1)
    K(i) = Parameter.m / 2 * SFrenet.Funb(SOL.Z(i,1))^2 * SOL.Z(i,2)^2;
    XYZ  = SFrenet.FunR (SOL.Z(i,1));
    V(i) = Parameter.m * Parameter.g * XYZ(3,1);
    E(i) = K(i) + V(i);

    FAux  = Parameter.m*Parameter.g * dot (E3,SFrenet.FunN(SOL.Z(i,1)));
    FN(i) = Parameter.m * SOL.Z(i,2)^2 * SFrenet.Funb(SOL.Z(i,1))^2 * SFrenet.FunKappa(SOL.Z(i,1)) + FAux;
    FAux  = Parameter.m*Parameter.g * dot (E3,SFrenet.FunB(SOL.Z(i,1)));
    FB(i) = FAux;
end

SOL.KE = K;
SOL.PE = V;
SOL.ME = E;
SOL.FN = FN;
SOL.FB = FB;

end

%%

function AUX = EtaDerivatives (SFrenet)

AUX.DDb    = simplify(diff (SFrenet.Db));
AUX.Dinvb  = simplify(diff(1.0/SFrenet.b));
AUX.DDx3   = SFrenet.DDR(3);
AUX.Dinvb2 = simplify(diff(1.0/SFrenet.b^2));

AUX.FunDDb    = matlabFunction(AUX.DDb);
if (nargin(AUX.FunDDb)==0)
    AUX.FunDDb = @(t) AUX.FunDDb();
end
AUX.FunDinvb  = matlabFunction(AUX.Dinvb);
if (nargin(AUX.FunDinvb)==0)
    AUX.FunDinvb = @(t) AUX.FunDinvb();
end
AUX.FunDDx3   = matlabFunction(AUX.DDx3);
if (nargin(AUX.FunDDx3)==0)
    AUX.FunDDx3 = @(t) AUX.FunDDx3();
end
AUX.FunDinvb2 = matlabFunction(AUX.Dinvb2);
if (nargin(AUX.FunDinvb2)==0)
    AUX.FunDinvb2 = @(t) AUX.FunDinvb2();
end

end