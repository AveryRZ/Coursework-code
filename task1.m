% 4-bus NR Power Flow (No SOP) with hand-written Jacobian
clear; clc;

% Slack bus
V1 = 1.01;             % |V1|, angle=0

% Specified injections (pu), generation +, load -
P2s = +0.85;  Q2s = 0.0;
P3s = -0.01;  Q3s = 0.0;
P4s = -0.5;   Q4s = -0.1;

% Branch impedances (pu, resistive)
Z12 = 0.04 + 0j;
Z23 = 0.02 + 0j;
Z14 = 0.04 + 0j;

Y12 = 1/Z12; Y23 = 1/Z23; Y14 = 1/Z14;

% Conductance coefficients (from Ybus)
G21 = -25; G22 = 75; G23 = -50;
G32 = -50; G33 = 50;
G41 = -25; G44 = 25;

% Pcalc and Qcalc for buses 2,3,4
Pcalc = @(V2,V3,V4,t2,t3,t4) [...
    V2*V1*G21*cos(t2) + V2^2*G22 + V2*V3*G23*cos(t2 - t3);...
    V3*V2*G32*cos(t3 - t2) + V3^2*G33;...
    V4*V1*G41*cos(t4) + V4^2*G44 ];

Qcalc = @(V2,V3,V4,t2,t3,t4) [...
    V2*V1*G21*sin(t2) + V2*V3*G23*sin(t2 - t3);...
    V3*V2*G32*sin(t3 - t2);...
    V4*V1*G41*sin(t4) ];

% Mismatch vector F(x) = [Pcalc-Pspec; Qcalc-Qspec] (6x1 column)
f = @(x) [...
    Pcalc(x(4),x(5),x(6),x(1),x(2),x(3)) - [P2s; P3s; P4s];...
    Qcalc(x(4),x(5),x(6),x(1),x(2),x(3)) - [Q2s; Q3s; Q4s] ];

% Hand-written Jacobian J(x): 6x6
J = @(x) jacobian_hand(V1, G21,G22,G23,G32,G33,G41,G44, x);

function Jm = jacobian_hand(V1,G21,G22,G23,G32,G33,G41,G44, x)
    t2=x(1); t3=x(2); t4=x(3); V2=x(4); V3=x(5); V4=x(6);

    % H = dP/dtheta
    H22 = 25*V1*V2*sin(t2) + 50*V2*V3*sin(t2 - t3);
    H23 = -50*V2*V3*sin(t2 - t3);
    H24 = 0;

    H32 = -50*V2*V3*sin(t3 - t2);
    H33 =  50*V2*V3*sin(t3 - t2);
    H34 = 0;

    H42 = 0; H43 = 0; H44 = 25*V1*V4*sin(t4);

    % N = dP/d|V|
    N22 = 150*V2 - 25*V1*cos(t2) - 50*V3*cos(t2 - t3);
    N23 = -50*V2*cos(t2 - t3);
    N24 = 0;

    N32 = -50*V3*cos(t3 - t2);
    N33 = 100*V3 - 50*V2*cos(t3 - t2);
    N34 = 0;

    N42 = 0; N43 = 0; N44 = 50*V4 - 25*V1*cos(t4);

    % M = dQ/dtheta
    M22 = -25*V1*V2*cos(t2) - 50*V2*V3*cos(t2 - t3);
    M23 =  50*V2*V3*cos(t2 - t3);
    M24 =  0;

    M32 =  50*V2*V3*cos(t3 - t2);
    M33 = -50*V2*V3*cos(t3 - t2);
    M34 =  0;

    M42 = 0; M43 = 0; M44 = -25*V1*V4*cos(t4);

    % L = dQ/d|V|
    L22 = -25*V1*sin(t2) - 50*V3*sin(t2 - t3);
    L23 = -50*V2*sin(t2 - t3);
    L24 = 0;

    L32 = -50*V3*sin(t3 - t2);
    L33 = -50*V2*sin(t3 - t2);
    L34 = 0;

    L42 = 0; L43 = 0; L44 = -25*V1*sin(t4);

    Jm = [...
        H22, H23, H24, N22, N23, N24;
        H32, H33, H34, N32, N33, N34;
        H42, H43, H44, N42, N43, N44;
        M22, M23, M24, L22, L23, L24;
        M32, M33, M34, L32, L33, L34;
        M42, M43, M44, L42, L43, L44 ];
end

% Newton iterations: x = [t2; t3; t4; V2; V3; V4]
x   = [0; 0; 0; 1.0; 1.0; 1.0];    % column vector
tol = 1e-10; maxIt = 100;

for k = 1:maxIt
    F  = f(x);                    % 6x1 column
    if norm(F, inf) < tol, break; end
    Jm = J(x);                    % 6x6
    dx = - Jm \ F;                % Newton step
    x  = x + dx;                  % update (column arithmetic)
    fprintf('Iter %2d, ||F||_inf = %.3e\n', k, norm(F, inf));
end

% Compose voltages
t2=x(1); t3=x(2); t4=x(3); V2=x(4); V3=x(5); V4=x(6);
V = [V1*exp(1j*0); V2*exp(1j*t2); V3*exp(1j*t3); V4*exp(1j*t4)];

% Branch currents
I12 = Y12*(V(1)-V(2));
I23 = Y23*(V(2)-V(3));
I14 = Y14*(V(1)-V(4));

fprintf('\n--- Solution (No SOP, hand-written Jacobian) ---\n');
for i=1:4
    fprintf('Bus %d: |V|=%.6f pu, angle=%.4f deg\n', i, abs(V(i)), rad2deg(angle(V(i))));
end
fprintf('Line currents: |I12|=%.5f, |I23|=%.5f, |I14|=%.5f pu\n', abs(I12), abs(I23), abs(I14));

inRange   = (abs(V)>=0.96) & (abs(V)<=1.04);
thermalOK = [abs(I12)<1.2, abs(I23)<1.2, abs(I14)<1.2];
fprintf('Voltage in-range flags: '); fprintf('%d ', inRange); fprintf('\n');
fprintf('Thermal check (<1.2 pu): %d %d %d\n', thermalOK);