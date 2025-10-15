clc; clear;
addpath(fullfile(pwd, 'functions'));


%% Simulation Parameters

c = physconst('LightSpeed'); % speed of light
fc = 25e6; % center frequency
lambda = c / fc; % wave length
kc = 2*pi/lambda; % wave number
nElem = 16; % number of array elements
dElem = 0.55*lambda; % element spacing

K = 1000; % fidelity of patterns
azAngles = linspace(-pi/2, pi/2, K).'; % azimuth scan angles
elAngles = linspace(-pi/2, pi/2, K).'; % elevation scan angles
coneAngles = asin(sin(azAngles) .* cos(elAngles));



%% Construct Array Geometry

% define array cartesian coordinates
elementXPositions = dElem * ((0: nElem-1) - (nElem-1)/2).';
elementYPositions = zeros(nElem, 1);
elementZPositions = zeros(nElem, 1);


%% Fmincon Adaptive Beamforming

% compute array response matrix across all coneAngles
S = exp(-1j * kc * sin(coneAngles) .* elementXPositions);

% define steer and null cone angles
phiSteer = [0];
phiNull = [40]; % [27:33];

% form steering vectors for steer and null angles
svSteer = exp(-1j * kc * sin(deg2rad(phiSteer)) .* elementXPositions);
svNull = exp(-1j * kc * sin(deg2rad(phiNull)) .* elementXPositions);
svSteerReal = [real(svSteer); image(svSteer)];

% define null depth objective function
minfun = @(w) (nullDepth(w, svNull));

% define nonlinear constraint function (forcing constant envelope)
nonlcon = @(w) (complexEnvelope(w));

% define inequality constraints
Aineq = [-real(svSteer)' -imag(svSteer)'];
bineq = [1.0];

% define equality constraints
Aeq = [-imag(svSteer)' real(svSteer)'];
beq = [0];

% define lower and upper filter weight bounds
lb = [];
ub = [];

% call fmincon (minimum of constrainted nonlinear multivariable function)
wopt = fmincon(minfun, svSteerReal, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
wopt = wopt(1:nElem) + 1j*wopt(nElem+1:2*nElem);

%  compute array factor output
arrayFactor = S' * wopt;
arrayFactor = arrayFactor / max(arrayFactor);

% plot array factor output



%% Functions

function [c, ceq] = complexEnvelope(w)
    c = [];
    n = length(w);
    x = w(1:n/2);
    y = w(n/2+1:n);
    ceq = x.^2 + y.^2 - ones(size(x));
end

function m = nullDepth(w, v_Null);
    n = length(w);
    x = w(1:n/2);
    y = w(n/2+1:n);
    m = max(abs((x + 1j*y)'*v_Null).^2);
end




