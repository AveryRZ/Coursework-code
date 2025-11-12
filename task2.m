% ============================================================
% SOP range via bisection: 4-bus NR with hand-written Jacobian
% Pt>0: draw at Bus 3, inject at Bus 4 (QSOP=0 both sides)
% ============================================================
clear; clc;

%% ----------- Network data (per-unit) -----------------------
V1  = 1.01;              % Slack magnitude (angle = 0 rad)
Z12 = 0.04 + 0j;
Z23 = 0.02 + 0j;
Z14 = 0.04 + 0j;

% Build Ybus directly
Ybus = zeros(4,4);
Y12  = 1/Z12;  Y23 = 1/Z23;  Y14 = 1/Z14;

% Branch 1-2
Ybus(1,1)=Ybus(1,1)+Y12;  Ybus(2,2)=Ybus(2,2)+Y12;
Ybus(1,2)=Ybus(1,2)-Y12;  Ybus(2,1)=Ybus(2,1)-Y12;
% Branch 2-3
Ybus(2,2)=Ybus(2,2)+Y23;  Ybus(3,3)=Ybus(3,3)+Y23;
Ybus(2,3)=Ybus(2,3)-Y23;  Ybus(3,2)=Ybus(3,2)-Y23;
% Branch 1-4
Ybus(1,1)=Ybus(1,1)+Y14;  Ybus(4,4)=Ybus(4,4)+Y14;
Ybus(1,4)=Ybus(1,4)-Y14;  Ybus(4,1)=Ybus(4,1)-Y14;

G = real(Ybus);   % B≈0 for this coursework
% Conductance entries used by hand-written formulas
G21 = G(2,1); G22 = G(2,2); G23 = G(2,3);
G32 = G(3,2); G33 = G(3,3);
G41 = G(4,1); G44 = G(4,4);

% Branch admittances for optional current checks
Y12 = 1/Z12; Y23 = 1/Z23; Y14 = 1/Z14;

%% ----------- Base injections (Monday 11am) -----------------
% Generation +, load -
P2s = +0.85;  Q2s = 0.0;
P3s0= -0.01;  Q3s = 0.0;
P4s0= -0.5;   Q4s = -0.1;

%% ----------- Define scalar |V3|(Pt) and endpoint values ----
V3_of_handle = @(Pt) V3_of(Pt, V1,G21,G22,G23,G32,G33,G41,G44, P2s,P3s0,P4s0,Q2s,Q3s,Q4s);

V3_minPt = V3_of_handle(-1);   % |V3| at Pt = -1
V3_maxPt = V3_of_handle(+1);   % |V3| at Pt = +1
fprintf('Endpoint |V3|: at Pt=-1 -> %.5f, at Pt=+1 -> %.5f\n', V3_minPt, V3_maxPt);

%% ----------- Bisection for upper limit (|V3|=1.04) ---------
% f_hi(Pt) = |V3(Pt)| - 1.04  (scalar function)
f_hi = @(Pt) V3_of_handle(Pt) - 1.04;
a_hi = -1; b_hi = +1;
fa = f_hi(a_hi); fb = f_hi(b_hi);

Pt_hi = NaN;
if fa*fb <= 0
    Pt_hi = bisect(f_hi, a_hi, b_hi, 1e-6, 60);
else
    if fa <= 0 && fb <= 0
        Pt_hi = -1;   % upper limit satisfied across entire range
    else
        warning('No crossing of 1.04 within [-1,1]; overvoltage persists for all Pt.');
    end
end

%% ----------- Bisection for lower limit (|V3|=0.96) ---------
% f_lo(Pt) = |V3(Pt)| - 0.96  (scalar function)
f_lo = @(Pt) V3_of_handle(Pt) - 0.96;
a_lo = -1; b_lo = +1;
fa = f_lo(a_lo); fb = f_lo(b_lo);

Pt_lo = NaN;
if fa*fb <= 0
    Pt_lo = bisect(f_lo, a_lo, b_lo, 1e-6, 60);
else
    % If V3(+1) >= 0.96, lower limit not hit within rating, so upper bound is Pt=+1
    if fb >= 0
        Pt_lo = +1;
    else
        % If V3(-1) <= 0.96, then even Pt=-1 already below lower limit
        Pt_lo = -1;
    end
end

%% ----------- Report the valid Pt range ---------------------
fprintf('\nValid Pt range (keep |V3| within [0.96,1.04]): [%.6f, %.6f] pu\n', Pt_hi, Pt_lo);

%% ----------- Optional: verify at the boundaries ------------
if ~isnan(Pt_hi)
    V_hi = NR_hand(V1,G21,G22,G23,G32,G33,G41,G44, P2s, P3s0-Pt_hi, P4s0+Pt_hi, Q2s, Q3s, Q4s);
    fprintf('At Pt=Pt_hi=%.6f: |V3|=%.5f\n', Pt_hi, abs(V_hi(3)));
end
if ~isnan(Pt_lo)
    V_lo = NR_hand(V1,G21,G22,G23,G32,G33,G41,G44, P2s, P3s0-Pt_lo, P4s0+Pt_lo, Q2s, Q3s, Q4s);
    fprintf('At Pt=Pt_lo=%.6f: |V3|=%.5f\n', Pt_lo, abs(V_lo(3)));
end

%% ----------- (Optional) show currents at a chosen Pt -------
Pt_inspect = 0.40;
V_chk = NR_hand(V1,G21,G22,G23,G32,G33,G41,G44, P2s, P3s0-Pt_inspect, P4s0+Pt_inspect, Q2s, Q3s, Q4s);
I12 = Y12*(V_chk(1)-V_chk(2));
I23 = Y23*(V_chk(2)-V_chk(3));
I14 = Y14*(V_chk(1)-V_chk(4));
fprintf('\nAt Pt=%.2f: |V|=[%.4f %.4f %.4f %.4f], |I12|=%.4f, |I23|=%.4f, |I14|=%.4f\n',...
    Pt_inspect, abs(V_chk(1)), abs(V_chk(2)), abs(V_chk(3)), abs(V_chk(4)), abs(I12), abs(I23), abs(I14));

% ===================== Local functions ======================
function v3 = V3_of(Pt, V1,G21,G22,G23,G32,G33,G41,G44, P2s,P3s0,P4s0,Q2s,Q3s,Q4s)
    % Return scalar |V3|(Pt)
    V = NR_hand(V1,G21,G22,G23,G32,G33,G41,G44,...
                P2s, P3s0 - Pt, P4s0 + Pt, Q2s, Q3s, Q4s);
    v3 = abs(V(3));
end

function Pt = bisect(fun, a, b, tol, itmax)
    fa = fun(a); fb = fun(b);
    if ~isscalar(fa) || ~isscalar(fb)
        error('fun(a) and fun(b) must be scalar values.');
    end
    if fa*fb > 0
        error('Bisection requires a sign change on [a,b].');
    end
    for it = 1:itmax
        c = 0.5*(a+b); fc = fun(c);
        if abs(fc) < 1e-10 || 0.5*(b-a) < tol
            Pt = c; return;
        end
        if fa*fc <= 0
            b = c; fb = fc;
        else
            a = c; fa = fc;
        end
    end
    Pt = 0.5*(a+b);
end

function V = NR_hand(V1,G21,G22,G23,G32,G33,G41,G44, P2s,P3s,P4s,Q2s,Q3s,Q4s)
    % Hand-written P,Q (B≈0). Variable x=[t2;t3;t4;V2;V3;V4]
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

    x   = [0;0;0; 1;1;1]; tol = 1e-10; maxIt = 100;
    for k=1:maxIt
        F=f(x); if norm(F,inf)<tol, break; end
        dx = - J(x) \ F; x = x + dx;
    end
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
    % Assemble Jacobian
    Jm = [...
        H22, H23, H24, N22, N23, N24;
        H32, H33, H34, N32, N33, N34;
        H42, H43, H44, N42, N43, N44;
        M22, M23, M24, L22, L23, L24;
        M32, M33, M34, L32, L33, L34;
        M42, M43, M44, L42, L43, L44 ];
end