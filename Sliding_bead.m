% Code for simulating the dynamic behavior of a sliding bead along an
% arbitrary curve in space. 
%--------------------------------------------------------------------------
% Author: Bruno Antonio Roccia
% Date: August 30, 2022
% Place: Bergen Offshore Wind Centre (BOW), Geophysical Institute, 
% University of Bergen, Norway
%--------------------------------------------------------------------------

% INPUT DATA

%       Parameter.m  : Mass
%       Parameter.g  : Gravitational acceleration
%       Parameter.Z1 : Initial condition (position)
%       Parameter.Z2 : Initial condition (velocity)
%       G1           : X-component of the curve parameterization
%       G2           : Y-component of the curve parameterization
%       G3           : Z-component of the curve parameterization

%       SOLVER

%       Solver.DT    : Time step
%       Solver.TOL   : Tolerance error
%       Solver.SimT  : Simulation time
%       Solver.Iter  : maximun number of iterations (implicit schemes)
%       Solver.IFLAG : Solver type --> ForwardEuler
%                                      BaackwardEuler
%                                      RungeKutta2
%                                      RungeKutta3
%                                      RungeKutta4
%                                      AdamsBashforth2
%                                      ODE15sMatlab
%                                      ODE23sMatlab, ODE23tMatlab, ODE23Matlab
%                                      ODE45Matlab
%                                      MidPointRule
%       Solver.FName : Name of the user-defined function 
%                      containing the Equations of Motion

%       PLOTS
%       Plot.IFLAG   : plot type --> Animation (motion visualization)
%                                    Energy
%                                    Portrait  (portrait diagram)
%                                    Force     (constraint forces)

%%

clear all
close all
clc

syms t real;

Parameter.m  = 1.0;
Parameter.g  = 9.81;
Parameter.Z1 = 2.6;
Parameter.Z2 = 0.0;

G1 = cos(t);
G2 = t^2*2;
G3 = sin(t);

% G1 = 2*cos(t)-1/2;
% G2 = sin(t);
% G3 = 3*cos(t)^2 - 2*cos(t) + 5/4;

Solver.DT     = 0.001;
Solver.TOL    = 1e-8;
Solver.SimT   = 5;
Solver.Iter   = 50;
Solver.IFLAG  = 'ForwardEuler'; 
Solver.FName  = 'Mysliding';

Plot.IFLAG    = 'Animation';
Plot.Interval = 15;
Plot.CInt   = 0.01;
Plot.CFinal = 4*pi;

%% PROBLEM 

R  = [G1; G2; G3];

[SFrenet] = SerretFrenet (R);

Parameter.Time = 0:Solver.DT:Solver.SimT;
Parameter.CI   = [Parameter.Z1; Parameter.Z2];

[SOL] = ODEIntegration (SFrenet,Solver,Parameter);

%% PLOT

if strcmpi(Plot.IFLAG, 'Portrait')
    plot (SOL.Z(:,1),SOL.Z(:,2),'color','b','linewidth',2)
    hold on
    h1 = plot (SOL.Z(1,1), SOL.Z(1,2),'o');
    set (h1,'markersize',8,'markerfacecolor','w','markeredgecolor','k','linewidth',2)
    grid on
    axis equal
    xlabel ('Z1')
    ylabel ('Z2')
    title ('Portrait diagram Z2 (velocity) vs Z1 (position)')
elseif strcmpi(Plot.IFLAG, 'Animation')
    CurveInterval = 0:Plot.CInt:Plot.CFinal;
    for i = 1:length(CurveInterval)
        XYZ(1:3,i) = SFrenet.FunR (CurveInterval(i));
    end
    plot3 (XYZ(1,:),XYZ(2,:),XYZ(3,:), 'color','b','linewidth',2);
    hold on
    plot3 (xlim, [0 0], [0 0],'color','k')
    plot3 ([0 0], ylim, [0 0],'color','k')
    plot3 ([0 0], [0 0], zlim,'color','k')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
    P = SFrenet.FunR (SOL.Z(1,1));
    h1 = plot3 (P(1), P(2), P(3), 'o');

   set (h1,'markersize',12,'markerfacecolor','g','markeredgecolor','k','linewidth',2)

    for i = 2:Plot.Interval:size(SOL.Z,1)
        P = SFrenet.FunR (SOL.Z(i,1));
        set (h1, 'xdata',P(1),'ydata',P(2),'zdata',P(3))
        plot3 (P(1),P(2),P(3),'marker','o','markersize',4,'markerfacecolor','r')
        pause(0.1)
    end
elseif strcmpi(Plot.IFLAG, 'Energy')
    plot (SOL.Tout,SOL.KE,'color','b','linewidth',2)
    hold on
    plot (SOL.Tout,SOL.PE,'color','r','linewidth',2)
    plot (SOL.Tout,SOL.ME,'color','k','linewidth',2)
    grid on
    xlabel ('Time [s]')
    ylabel ('Energy')
    legend ('Kinetic Energy', 'Potential Energy', 'Mechanical Energy')
    title('Energy of the system')
elseif strcmpi(Plot.IFLAG, 'Force')
    plot (SOL.Tout,SOL.FN,'color','b','linewidth',2)
    hold on
    plot (SOL.Tout,SOL.FB,'color','r','linewidth',2)
    grid on
    xlabel ('Time [s]')
    ylabel ('Constraint Forces')
    legend ('F along N-direction', 'F along B direction')
    title('Constraint forces along the normal and binormal directions')
elseif strcmpi(Plot.IFLAG, 'TimeSeries')
    plot (SOL.Tout,SOL.Z(:,1),'color','b','linewidth',2)
    hold on
    plot (SOL.Tout,SOL.Z(:,2),'color','r','linewidth',2)
    grid on
    xlabel ('Time [s]')
    ylabel ('Z1 and Z2')
    legend ('Position Z1', 'Velocity Z2')
    title('Time series response')
end

%%

function [SFrenet] = SerretFrenet (R)

DR    = simplify(diff (R));
b     = sqrt (dot(DR,DR));

T     = simplify(DR / b);
Db    = diff(b);

DDR   = simplify(diff (DR));
PT    = eye(3) - T*transpose(T);

AUX   = norm (cross(DDR,T),2);

N     = PT*DDR / AUX;
B     = 1/AUX * cross(T,DDR);

MDT   = simplify(1/b * AUX);

Kappa = simplify(1/b^3 * AUX);

SFrenet.R     = R;
SFrenet.DR    = DR;
SFrenet.DDR   = DDR;
SFrenet.b     = b;
SFrenet.Db    = Db;
SFrenet.T     = T;
SFrenet.N     = N;
SFrenet.B     = B;
SFrenet.MDT   = MDT;
SFrenet.Kappa = Kappa;

SFrenet.FunR     = matlabFunction(R);
SFrenet.FunDR    = matlabFunction(DR);
SFrenet.Funb     = matlabFunction(b);
if (nargin(SFrenet.Funb)==0)
    SFrenet.Funb = @(t) SFrenet.Funb();
end
SFrenet.FunDb = matlabFunction(Db);
if (nargin(SFrenet.FunDb)==0)
    SFrenet.FunDb = @(t) SFrenet.FunDb();
end
SFrenet.FunT = matlabFunction(T);
if (nargin(SFrenet.FunT)==0)
    SFrenet.FunT = @(t) SFrenet.FunT();
end
SFrenet.FunN = matlabFunction(N);
if (nargin(SFrenet.FunN)==0)
    SFrenet.FunN = @(t) SFrenet.FunN();
end
SFrenet.FunB = matlabFunction(B);
if (nargin(SFrenet.FunB)==0)
    SFrenet.FunB = @(t) SFrenet.FunB();
end
SFrenet.FunMDT = matlabFunction(MDT);
if (nargin(SFrenet.FunMDT)==0)
    SFrenet.FunMDT = @(t) SFrenet.FunMDT();
end
SFrenet.FunKappa = matlabFunction(Kappa);
if (nargin(SFrenet.FunKappa)==0)
    SFrenet.FunKappa = @(t) SFrenet.FunKappa();
end

end
