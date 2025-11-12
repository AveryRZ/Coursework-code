% ============================================================
% Contingency: Pt=0.4 (draw at Bus 3, inject at Bus 4), Bus 2 generator trips
% 4-bus NR power flow with hand-written Jacobian
% Generation positive, load negative; Slack bus: Bus 1 (|V1|=1.01, angle=0)
% ============================================================
clear; clc;

%% ----------- Network data (per-unit, resistive) ------------
V1  = 1.01;              % Slack magnitude, angle = 0
Z12 = 0.04 + 0j;
Z23 = 0.02 + 0j;
Z14 = 0.04 + 0j;

% Build Ybus directly
Ybus = zeros(4,4);
Y12 = 1/Z12; Y23 = 1/Z23; Y14 = 1/Z14;

% Branch 1-2
Ybus(1,1)=Ybus(1,1)+Y12; Ybus(2,2)=Ybus(2,2)+Y12;
Ybus(1,2)=Ybus(1,2)-Y12; Ybus(2,1)=Ybus(2,1)-Y12;
% Branch 2-3
Ybus(2,2)=Ybus(2,2)+Y23; Ybus(3,3)=Ybus(3,3)+Y23;
Ybus(2,3)=Ybus(2,3)-Y23; Ybus(3,2)=Ybus(3,2)-Y23;
% Branch 1-4
Ybus(1,1)=Ybus(1,1)+Y14; Ybus(4,4)=Ybus(4,4)+Y14;
Ybus(1,4)=Ybus(1,4)-Y14; Ybus(4,1)=Ybus(4,1)-Y14;

G = real(Ybus);                % B≈0 for this coursework

% Conductance entries used by hand-written formulas
G21 = G(2,1); G22 = G(2,2); G23 = G(2,3);
G32 = G(3,2); G33 = G(3,3);
G41 = G(4,1); G44 = G(4,4);

%% ----------- Scenario: Pt=0.4 and generator at Bus 2 trips --
Pt  = 0.4;

% Specified injections (generation +, load -)
P2s = 0.0;   Q2s = 0.0;             % generator tripped
P3s = -0.01 - Pt;  Q3s = 0.0;       % SOP draws at Bus 3
P4s = -0.5  + Pt;  Q4s = -0.1;      % SOP feeds at Bus 4

%% ----------- Solve with hand-written NR ---------------------
[V, iter] = NR_hand(V1, G21,G22,G23,G32,G33,G41,G44, P2s,P3s,P4s,Q2s,Q3s,Q4s);

% Branch currents and checks
I12 = Y12*(V(1)-V(2));
I23 = Y23*(V(2)-V(3));
I14 = Y14*(V(1)-V(4));

fprintf('--- Contingency result (Pt=0.4, Bus2 trip) ---\n');
for i=1:4
    fprintf('Bus %d: |V|=%.6f pu, angle=%.4f deg\n', i, abs(V(i)), rad2deg(angle(V(i))));
end
fprintf('Line currents: |I12|=%.5f, |I23|=%.5f, |I14|=%.5f pu\n', abs(I12), abs(I23), abs(I14));

% Statutory voltage range check
inRange = (abs(V)>=0.96) & (abs(V)<=1.04);
fprintf('Voltage statutory in-range flags [B1 B2 B3 B4]: ');
fprintf('%d ', inRange); fprintf('\n');

% Thermal check (<1.2 pu)
thermalOK = [abs(I12)<1.2, abs(I23)<1.2, abs(I14)<1.2];
fprintf('Thermal check (<1.2 pu): %d %d %d\n', thermalOK);

% ===================== Local functions ======================
function [V, iter] = NR_hand(V1,G21,G22,G23,G32,G33,G41,G44,...
                             P2s,P3s,P4s,Q2s,Q3s,Q4s)
    % Hand-written P,Q (B≈0)
    Pcalc = @(V2,V3,V4,t2,t3,t4) [...
        V2*V1*G21*cos(t2) + V2^2*G22 + V2*V3*G23*cos(t2 - t3);...
        V3*V2*G32*cos(t3 - t2) + V3^2*G33;...
        V4*V1*G41*cos(t4) + V4^2*G44 ];
    Qcalc = @(V2,V3,V4,t2,t3,t4) [...
        V2*V1*G21*sin(t2) + V2*V3*G23*sin(t2 - t3);...
        V3*V2*G32*sin(t3 - t2);...
        V4*V1*G41*sin(t4) ];

    f = @(x) [...
        Pcalc(x(4),x(5),x(6),x(1),x(2),x(3)) - [P2s; P3s; P4s];...
        Qcalc(x(4),x(5),x(6),x(1),x(2),x(3)) - [Q2s; Q3s; Q4s] ];

    J = @(x) jacobian_hand(V1, G21,G22,G23,G32,G33,G41,G44, x);

    % Newton iterations
    x   = [0; 0; 0; 1.0; 1.0; 1.0];  % flat start
    tol = 1e-10; maxIt = 100;
    for iter = 1:maxIt
        F  = f(x);
        if norm(F, inf) < tol, break; end
        Jm = J(x);
        dx = - Jm \ F;
        x  = x + dx;
    end

    % Compose complex voltages
    t2=x(1); t3=x(2); t4=x(3); V2=x(4); V3=x(5); V4=x(6);
    V = [V1*exp(1j*0); V2*exp(1j*t2); V3*exp(1j*t3); V4*exp(1j*t4)];
end

function Jm = jacobian_hand(V1,G21,G22,G23,G32,G33,G41,G44, x)
    t2=x(1); t3=x(2); t4=x(3); V2=x(4); V3=x(5); V4=x(6);

    % H = dP/dtheta
    H22 = -V2*V1*G21*sin(t2) - V2*V3*G23*sin(t2 - t3);
    H23 =  +V2*V3*G23*sin(t2 - t3);
    H24 = 0;

    H32 =  +V2*V3*G32*sin(t3 - t2);
    H33 =  -V2*V3*G32*sin(t3 - t2);
    H34 = 0;

    H42 = 0; H43 = 0;
    H44 = -V4*V1*G41*sin(t4);

    % N = dP/d|V|
    N22 = 2*V2*G22 + V1*G21*cos(t2) + V3*G23*cos(t2 - t3);
    N23 = V2*G23*cos(t2 - t3);
    N24 = 0;

    N32 = V3*G32*cos(t3 - t2);
    N33 = 2*V3*G33 + V2*G32*cos(t3 - t2);
    N34 = 0;

    N42 = 0; N43 = 0;
    N44 = 2*V4*G44 + V1*G41*cos(t4);

    % M = dQ/dtheta
    M22 =  V2*V1*G21*cos(t2) + V2*V3*G23*cos(t2 - t3);
    M23 = -V2*V3*G23*cos(t2 - t3);
    M24 = 0;

    M32 = -V2*V3*G32*cos(t3 - t2);
    M33 =  V2*V3*G32*cos(t3 - t2);
    M34 = 0;

    M42 = 0; M43 = 0;
    M44 =  V4*V1*G41*cos(t4);

    % L = dQ/d|V|
    L22 =  V1*G21*sin(t2) + V3*G23*sin(t2 - t3);
    L23 =  V2*G23*sin(t2 - t3);
    L24 = 0;

    L32 =  V3*G32*sin(t3 - t2);
    L33 =  V2*G32*sin(t3 - t2);
    L34 = 0;

    L42 = 0; L43 = 0;
    L44 =  V1*G41*sin(t4);

    % Assemble Jacobian: rows [P2 P3 P4 Q2 Q3 Q4], cols [t2 t3 t4 V2 V3 V4]
    Jm = [...
        H22, H23, H24, N22, N23, N24;
        H32, H33, H34, N32, N33, N34;
        H42, H43, H44, N42, N43, N44;
        M22, M23, M24, L22, L23, L24;
        M32, M33, M34, L32, L33, L34;
        M42, M43, M44, L42, L43, L44 ];
end