function [DZ] = Mysliding(t,Z, Parameter, SFrenet)

DZ = zeros(2,1);

Db  = SFrenet.FunDb(Z(1,1));
b   = SFrenet.Funb(Z(1,1));
DX3 = SFrenet.FunDR(Z(1,1));
DX3 = DX3(3,1);
g   = Parameter.g;

DZ(1,1) = Z(2,1);
DZ(2,1) = -Db/b * Z(2,1)^2 - g*DX3/b^2;

end