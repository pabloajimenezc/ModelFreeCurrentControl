function cost = simulate(x)

rng(0)

% plecs
% open('Simulate.slx');

% System parameters
% Frequencies
fs = 40e3;   % Sampling frequency
Ts = 1/fs;   % Sampling step
Ti = Ts/10; % Simulation step
fc = fs;     % Carrier frequency
f = 50;      % Fundamental output frequency
w = 2*pi*f;
% Converter
Vdc = 700;
Cul = 2200e-6;
% Filter
Rf = 1e-3;
Lf = 2e-3;
Cf = 50e-6;
Vo = 230;
% Load
Ro = 180;
Lo = 2e-3;
% Measurements
% Current
i_max = 50;             % Saturation
i_std = 0.01 * i_max / 3; % White noise with 3-sigma rule
% Voltage
v_max = 1000;            % Saturation
v_std = 0.01 * v_max / 3; % White noise with 3-sigma rule

% System model
% Equilibrium equations

% Continuous time


% Model matrices
I2 = eye(2);
O2 = zeros(2);

A = [-Rf/Lf*I2, -1/Lf*I2;
       1/Cf*I2,       O2];
B = [Vdc/(3*Lf)*I2;
                O2];
E = [      O2;
     -1/Cf*I2];

% Discrete time
% For control
% One step

nx = size(A, 2); % # of state variables
nu = size(B, 2); % # of control actions
np = size(E, 2); % # of disturbances

Ad = expm(Ts*A);
integrated = A\(Ad-eye(nx));
Bd = integrated * B;
Ed = integrated * E;

% Np steps
% For control
Np = 5;

% With

AA = arrayfun(@(t) Ad^t, 1:Np, 'UniformOutput', false);
AA = vertcat(AA{:});

% AAA
AA_aux = [eye(nx); AA(1:end-nx, :)];
AAA = cell2mat(arrayfun(@(t) [zeros(nx*t, nx); AA_aux(1:end-nx*t, :)], 0:Np-1, 'UniformOutput', false));
% IB
IB = repmat({Bd}, 1, Np);
IB = blkdiag(IB{:});
% BB
BB = AAA * IB;
% IE
IE = repmat({Ed}, 1, Np);
IE = blkdiag(IE{:});
% EE
EE = AAA * IE;

% Kalman Filter
% For current load estimation
% Negative and positive sequence, alpha-beta components for each harmonic
% That is 4 estimated variables for each harmonic
nh = 3;                                      % # of estimated harmonics
H_vec = find(mod(1:100, 2) & mod(1:100, 3)); % Odds but non multiples of 3
H_vec = H_vec(1:nh);

nH = 4 * nh;   % # of KF state variables
J     = [0, -1;  % 90° rotation matrix
         1,  0];
JJ    = [J,  O2; % 90° and -90° rotation matrix
         O2, -J];
% State matrix
A_KF = arrayfun(@(h) h*w*JJ, H_vec, 'UniformOutput', false);
A_KF = blkdiag(A_KF{:});
Ad_KF = expm(Ts * A_KF);
% Output matrix
C_KF = [1, 0, 1, 0;
        0, 1, 0, 1];
C_KF = repmat(C_KF, 1, nh);
% Kalman gain
r = i_std^2;
R_KF = r * I2;
q = r * x.q;
Q_KF = q * eye(nH);
[P_KF, KF, L_KF] = idare(Ad_KF.', C_KF.', Q_KF, R_KF);
KF = KF.';

% Control
% Grid-forming capacitor voltage control
% Cost function

% Reference tracking
% For  and 
% Weights
lambda_i = x.lambda_i; % Current
Wx = eye(nx);
Wx(1, 1) = lambda_i / i_max^2;
Wx(2, 2) = lambda_i / i_max^2;
Wx(3, 3) = 1 / v_max^2;
Wx(4, 4) = 1 / v_max^2;

Wx = repmat({Wx}, 1, Np);
Wx = blkdiag(Wx{:});
% Cost function

% With

% Then

% Finally,

Hx = 2 * (BB.') * Wx * BB;
% fx = (Wx * BB).' * (-X_ref + AA * x0 + EE * P);

% Control effort
% Weights
lambda_u = x.lambda_u;
% Cost function

% Finally,

Hu = 2 * eye(np*Np) / v_max^2;
% fu = -2 * U_ref / v_max^2;

% Total cost function
H = Hx + lambda_u * Hu;
H = (H + H.') / 2;
% f = fx + lambda_u * fu;

% Steady-state references
% Output voltage

% Converter current from the LCK
% Converter voltage from the LVK
% From the system equations
% Converter current
% 
% Converter voltage
% 
% In vectorial form
% 
% 
% Control action constraints
% Control actions
% Active phasor hexagon


Aineq = [sqrt(3)   1
         0         1
        -sqrt(3)   1
        -sqrt(3)  -1
         0        -1
         sqrt(3)  -1];

bineq = sqrt(3) * [1
                   1/2
                   1
                   1
                   1/2
                   1];

AAineq = repmat({Aineq}, 1, Np);
AAineq = blkdiag(AAineq{:});

bbineq = repmat({bineq}, 1, Np);
bbineq = vertcat(bbineq{:});

% Predictions
% Assuming steady state for every prediction
% Also, every vector of predictions has the actual variable value as first element
rot = @(g) [cos(g), -sin(g); sin(g), cos(g)]; % Rotation matrix
% Control action reference
u_ref_ROT = cell2mat(arrayfun(@(t) rot(w*Ts*t), (0:Np-1)', 'UniformOutput', false));
% State variables references
x_ref_ROT = cell2mat(arrayfun(@(t) [rot(w*Ts*t), O2; O2, rot(w*Ts*t)], (0:Np-1)', 'UniformOutput', false));
% Disturbance (load current)

% Every harmonic  in  coordinates is predicted rotating it an angle , counter clockwise for the positive sequence and clockwise for the negative sequence.

pn_rot = @(h, t) [rot(H_vec(h)*w*Ts*t), O2; O2, rot(-H_vec(h)*w*Ts*t)]; % Positive and negative sequence rotation
Ht_rot = @(t) arrayfun(@(h) pn_rot(h, t), (1:nh)', 'UniformOutput', false);
HT_ROT = cell(Np, 1);
for t = 1:Np
    Ht = Ht_rot(t-1);
    HT_ROT{t} = blkdiag(Ht{:});
end
HT_ROT = cell2mat(HT_ROT);
x_KF_ROT = HT_ROT;
% Then to obtain y_KF

CC_KF = repmat({C_KF}, 1, Np);
CC_KF = blkdiag(CC_KF{:});

KF_ROT = CC_KF * x_KF_ROT;
dKF_ROT = CC_KF * x_KF_ROT * A_KF;

% Modulation
% M = 
% 1: DPWMmin
% 2: DPWmax
% 3: SVPWM
% 4: DPWM0
% 5: DPWM1
% 6: DPWM2
% 7: DPWM3
% 8: SPWM TAMALA
% HB: Capacitor voltage balancing band

% Simulation
Tf = 3 / f; %[s] Total simulation time
% out = sim('Simulate.slx');

vars = who;
for k = 1:length(vars)
    assignin('base', vars{k}, eval(vars{k}));
end
out = sim('Simulate');

cost = 1;

end