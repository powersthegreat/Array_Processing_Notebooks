clc; clear;


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



%% Construct LPDA Element Geometry

% % debugging for faster simulations
% element = design(dipole, fc);

% configure lpda structure geometry
element = lpda;
element.BoardLength = 0.75 * lambda;
element.BoardWidth = 0.75 * lambda;
element.Height = 0.01 * lambda;
element.StripLineWidth = 0.005 * lambda;
element.FeedLength = 0.025 * lambda;

% configure lpda arm geometry
nArms = 8;
armStart = 0.16 * lambda;
armEnd = 0.34 * lambda;
element.ArmLength = linspace(armStart, armEnd, nArms);
element.ArmWidth = 0.02 * element.ArmLength;
element.ArmSpacing = linspace(0.08, 0.12, nArms-1) * lambda;

% rotate elements along y axis
element.TiltAxis = [0 1 0; 1 1 0];
element.Tilt = [180 180];

% plot single element geometry
figure("Name", "Single Element Geometry");
show(element);
title("LPDA Antenna Configured for 25MHz");
axis equal;
grid on;
view(0, 270);

% % plot single element impedance (HIGH COMPUTATION!)
% figure("Name", "Single Element Impedance Sweep");
% freqSweep = linspace(fc - (fc*.1), fc + (fc*.1), 10);
% impedance(element, freqSweep);
% 
% % plot single element bandwidth (HIGH COMPUTATION!)
% figure("Name", "Single Element Bandwidth");
% freqSweep = linspace(fc - (fc*.1), fc + (fc*.1), 10);
% absBW = bandwidth(element, freqSweep);



%% Compute Single Element Patterns (Az and El)

elementPatternAz = 10.^(pattern(element, fc, rad2deg(azAngles-pi/2), 0, 'Type', 'gain')/10).';
elementPatternEl = 10.^(pattern(element, fc, 0, rad2deg(azAngles-pi/2), 'Type', 'gain')/10).';
load('data/elementPatternAz.mat');
load('data/elementPatternEl.mat');

% normalize element pattern
elementPatternAz = elementPatternAx / max(elementPatternAz);
elementPatternEl = elementPatternEl / max(elementPatternEl);

% save element patterns for later
save('data/elementPatternAz.mat', 'elementPatternAz');
save('data/elementPatternEl.mat', 'elementPatternEl');

% plot element azimuth and elevation patterns
figure("Name", "Single Element Azimuth and Elevation Patterns", "Color", "w");
subplot(2, 1, 1);
plot(rad2deg(azAngles), 10*log10(elementPatternAz.^2), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2.0);
xlabel('Azimuth Angle (deg)');
ylabel('Gain (dB)');
title('Element Azimuth Pattern');
grid on;
xlim([-90 90]);
subplot(2, 1, 2);
plot(rad2deg(azAngles), 10*log10(elementPatternEl.^2), 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 2.0);
xlabel('Elevation Angle (deg)');
ylabel('Gain (dB)');
title('Element Elevation Pattern');
grid on;
xlim([-90 90]);
sgtitle('Single Element Radiation Patterns (Normalized Gain)');
set(gcf, 'Position', [100, 100, 700, 600]);



%% Construct Array Geometry

% define array cartesian coordinates
elementXPositions = dElem * ((0: nElem-1) - (nElem-1)/2).';
elementYPositions = zeros(nElem, 1);
elementZPositions = zeros(nElem, 1);

% generate confromal array geometry
array = conformalArray;
array.Element = element;
array.ElementPosition = [elementXPositions, elementYPositions, elementZPositions];

% plot 3-dimensional array geometry
figure("Name", "Array Geometry");
show(array);



%% Compute Array Factors (Az and El)

% define steering angles
azSteer = deg2rad(0);
elSteer = deg2rad(0);

% convert steering angles to direction vector (cartesian)
u = [sin(azSteer)*cos(elSteer), cos(azSteer)*cos(elSteer), sin(elSteer)].';

% form steering vector (and normalize)
svSteer = exp(-1j * kc * (array.ElementPosition * u));
svSteer = svSteer / norm(svSteer);

% form array response matrix across fixed elevation
UFixedEl = [sin(azAngles)*cos(zeros(K, 1)), cos(azAngles)*cos(zeros(K, 1)), sin(zeros(K, 1))].';
S = exp(-1j * kc * (array.ElementPosition * UFixedEl));

% compute 2D array factor across azimuth sweep
arrayFactorAz = S' * svSteer;
arrayFactorAz = arrayFactorAz / maz(arrayFactorAz);

% form array response matrix across fixed azimuth
UFixedAz = [sin(zeros(K, 1))*cos(elAngles), cos(zeros(K, 1))*cos(elAngles), sin(elAngles)].';
S = exp(-1j * kc * (array.ElementPosition * UFixedAz));

% compute 2D array factor across elevation sweep
arrayFactorEl = S' * svSteer;
arrayFactorEl = arrayFactorEl / max(arrayFactorEl);

% plot 2D array factor across azimuth
figure("Name", "Array Factor - Azimuth Sweep", "Color", "w");
plot(rad2deg(azAngles), 20*log10(abs(arrayFactorAz)), 'Color', [0.13, 0.55, 0.13], 'LineWidth', 2.0);
xlabel('Azimuth Angle (deg)');
ylabel('Normalized Gain (dB)');
title('2D Array Factor vs. Azimuth');
grid on;
xlim([-90 90]);
ylim([-60 0]);
set(gca, 'FontSize', 12);


% % plot 2D array factor across elevation
% figure("Name", "Array Factor - Elevation Sweep", "Color", "w");
% plot(rad2deg(elAngles), 20*log10(abs(arrayFactorEl)), 'Color', [0.85, 0.33, 0.10], 'LineWidth', 2.0);
% xlabel('Elevation Angle (deg)');
% ylabel('Normalized Gain (dB)');
% title('2D Array Factor vs. Elevation');
% grid on;
% xlim([-90 90]);
% ylim([-60 0]);
% set(gca, 'FontSize', 12);

% % compute array factor over azimuth and elevation
% [AZ, EL] = meshgrid(azAngles, elAngles);
% Ux = sin(AZ) .* cos(EL);
% Uy = cos(AZ) .* cos(EL);
% Uz = sin(EL);
% U = [Ux(:); Uy(:); Uz(:)];
% S = exp(-1j * kc * (array.ElementPosition * U));
% AF = svSteer' * S;
% AF = reshape(AF, K, K);
% AF_dB = 20*log10(abs(AF) / max(abs(AF(:))));

% % plot #D array factor across both azimuth and elevation
% figure("Name", "3D Array Factor", "Color", "w");
% imagesc(rad2deg(azAngles), rad2deg(elAngles), AF_dB);
% set(gca, 'YDir', 'normal');
% xlabel('Azimuth Angle (deg)');
% ylabel('Elevation Angle (deg)');
% title('3D Array Factor (Azimuth vs Elevation)');
% colorbar;
% colormap(turbo);
% caxis([-40 0]);
% set(gca, 'FontSize', 12);



%% Compute Array Patterns (Az and El)

arrayPatternAz = elementPatternAz .* arrayFactorAz;
arrayPatternEl = elementPatternEl .* arrayFactorEl;
arrayPatternAz = arrayPatternAz / max(abs(arrayPatternAz));
arrayPatternEl = arrayPatternEl / max(abs(arrayPatternEl));

% plot array pattern across fixed elevation
figure("Name", "Array Pattern - Azimuth Sweep", "Color", "w");
plot(rad2deg(azAngles), 20*log10(abs(arrayPatternAz)), 'Color', [0.13, 0.55, 0.13], 'LineWidth', 2.0);
xlabel('Azimuth Angle (deg)');
ylabel('Normalized Gain (dB)');
title('Array Pattern vs. Azimuth');
grid on;
xlim([-90 90]);
ylim([-60 0]);
set(gca, 'FontSize', 12);

% plot array pattern across fixed aximuth
figure("Name", "Array Pattern - Elevation Sweep", "Color", "w");
plot(rad2deg(elAngles), 20*log10(abs(arrayPatternEl)), ...
    'Color', [0.85, 0.33, 0.10], 'LineWidth', 2.0);
xlabel('Elevation Angle (deg)');
ylabel('Normalized Gain (dB)');
title('Array Pattern vs. Elevation');
grid on;
xlim([-90 90]);
ylim([-60 0]);
set(gca, 'FontSize', 12);

% compute 3D pattern across both azimuth and elevation
[AZ, EL] = meshgrid(azAngles, elAngles);
Ux = sin(AZ) .* cos(EL);
Uy = cos(AZ) .* cos(EL);
Uz = sin(EL);
U = [Ux(:)'; Uy(:)'; Uz(:)'];
S = exp(-1j * kc * (array.ElementPosition * U));
AF = svSteer' * S;
AF = reshape(AF, K, K);

% create 2D element pattern (outer product)
EP = elementPatternAz(:) .* elementPatternEl(:).'; 

% compute total array pattern
AP = AF .* EP;
AP_dB = 20 * log10(abs(AP) / max(abs(AP(:))));

% plot 3D array pattern as heat map
figure("Name", "3D Array Pattern (Az vs El)", "Color", "w");
imagesc(rad2deg(azAngles), rad2deg(elAngles), AP_dB);
set(gca, 'YDir', 'normal');
xlabel('Azimuth Angle (deg)');
ylabel('Elevation Angle (deg)');
title('3D Array Pattern (Azimuth vs Elevation)');
colorbar;
colormap(turbo);
clim([-40 0]);
set(gca, 'FontSize', 12);

% plot 3D array pattern as mesh surface
figure("Name", "3D Surface Array Pattern", "Color", "w");
surf(rad2deg(AZ), rad2deg(EL), AP_dB, 'EdgeColor', 'none');
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
zlabel('Gain (dB)');
title('3D Array Pattern Surface');
colormap(turbo);
colorbar;
caxis([-40 0]);
view(45, 30);
grid on;
set(gca, 'FontSize', 12);



%% Simulate Mutual Coupling Between Elements

% compute s-parameters using matlabs built in solver
sObj = sparameters(arrayt, fc, mean(real(impedance(Array, fc))));

% plot s-parameter matrix
figure("Name", "Array Mutual Coupling Matrix");
imagesc(db(sObj.Parameters + eye(nElem)));

% compute active reflected power across elements
reflectionCeofs = zeros(nElem, K);

for idx = 1:K

    % define steering angle and form steering vector
    azSteer = azAngles(idx);
    u = [sin(azSteer)*cos(0), cos(azSteer)*cos(0), sin(0)].';
    svSteer = exp(-1j * kc * (array.ElementPosition * u));

    % form incident and reflected voltages at each element
    sIncident = svSteer;
    sReflected = sObj.Parameters * svSteer;

    % compute reflection coefficients and store
    reflectionCeofs(:, idx) = sReflected ./ sIncident;

end

% plot reflected power across elements
reflectionCoefsdB = db(abs(reflectionCeofs).^2);
figure("Name", "Active Reflected Power Across Array");
imagesc(elementXPositions/dElem, rad2deg(azAngles), reflectionCoefsdB.');

% compute reflected power against thresold (NO MORE THAN 3dB)
maxVSWR = 3;
maxReflect = db(abs((maxVSWR - 1) / (maxVSWR + 1)));
reflectionCoefsdB(reflectionCoefsdB > maxReflect) = -Inf;

% plot reflected power against threshold
figure("Name", "Active Reflected Power Across Array Against Threshold");
imagesc(elementXPositions/dElem, rad2deg(azAngles), reflectionCeofsdB.');



% %% Compute Active VSWR Across Elements

% VSWRCoeffs = zeros(nElem, K);

% for idx = 1:K

%     % define steering angle and form steering vector
%     azSteer = azAngles(idx);
%     u = [sin(azSteer)*cos(0), cos(azSteer)*cos(0), sin(0)].';
%     svSteer = exp(-1j * kc * (array.ElementPosition * u));

%     % form incident and reflected voltages at each element
%     sIncident = svSteer;
%     sReflected = sObj.Parameters * svSteer;

%     % compute VSWR and store
%     reflectionCoef = sReflected ./ sIncident;
%     VSWRCoeffs(:, idx) = (1 + abs(reflectionCeof)) ./ (1 - abs(reflectionCeof));

% end

% % plot VSWR sweep
% figure("Name", "VSWR Sweep");
% iamgesc(rad2deg(azAngles), elementXPositions/dElem, VSWRCoeffs);

% % plot usable VSWR across array
% figure("Name", "Threholded Active VSWR");
% clim([0 min(VSWRCoeffs(:)) + 10]);
% imagesc(rad2deg(azAngles), elementXPositions/dElem, VSWRCoeffs);



% %% Adaptive Transmit Beamforming w/ Coupling

% % define steering angle and form steering vector
% azSteer = azAngles(0);
% u = [sin(azSteer)*cos(0), cos(azSteer)*cos(0), sin(0)].';
% svSteer = exp(-1j * kc * (array.ElementPosition * u));

% % define null angle and form nulling vector
% azNull = azAngles(35);
% u = [sin(azNull)*cos(0), cos(azNull)*cos(0), sin(0)].';
% svNull = exp(-1j * kc * (array.ElementPosition * u));

% % form interference covariance matrix (and diagonally load)
% dlFactor = 10^(-30/20);
% R = (svNull * svNull') + (dlFactor .* eye(nElem));

% % form array response matrix across fixed elevation
% UFixedEl = [sin(azAngles)*cos(zeros(K, 1)), cos(azAngles)*cos(zeros(K, 1)), sin(zeros(K, 1))].';
% S = exp(-1j * kc * (array.ElementPosition * UFixedEl));

% % solve for MVDR filter coefficients
% gain = max(S' * svSteer);
% wMVDR = (gain * R^-1 * svSteer) / (svSteer' * R^-1 * svSteer);

% % compute 2D array factor across azimuth sweep
% arrayFactorAz = S' * wMVDR;
% arrayFactorCoupledAz = S' * (wMVDR .* ((eye(nElem) + sObj.Parameters) * ones(nElem, 1)));

% % form array response matrix across fixed azimuth
% UFixedAz = [sin(zeros(K, 1))*cos(elAngles), cos(zeros(K, 1))*cos(elAngles), sin(elAngles)].';
% S = exp(-1j * kc * (array.ElementPosition * UFixedAz));

% % compute 2D array factor across elevation sweep
% arrayFactorEl = S' * wMVDR;
% arrayFactorCoupledEl = S' * (wMVDR .* ((eye(nElem) + sObj.Parameters) * ones(nElem, 1)));

% % plot 2D array factor across azimuth

% % plot 2D array factor across elevation

% % compute and plot 3D array factor across azimuth and elevation



% %% Compute Array Patterns with Coupling



% %% Simulate Array Perturbatios

% % define phase noise steering vector
% sigma = 10^(-35/10);
% svPhaseNoise = exp(-1j * swrt(sigma) * randn(nElem, 1));

% % define gain noise steering vector
% sigma = 10^(-40/10);
% svGainNoise = ones(nElem, 1) + (swrt(signma) * randn(nElem, 1));



% %% Adaptive Beamforming with Perturbations

% % define steering angle and form steering vector
% azSteer = azAngles(0);
% u = [sin(azSteer)*cos(0), cos(azSteer)*cos(0), sin(0)].';
% svSteer = exp(-1j * kc * (array.ElementPosition * u));

% % define null angle and form nulling vector
% azNull = azAngles(35);
% u = [sin(azNull)*cos(0), cos(azNull)*cos(0), sin(0)].';
% svNull = exp(-1j * kc * (array.ElementPosition * u));

% % form interference covariance matrix (and diagonally load)
% dlFactor = 10^(-30/20);
% R = (svNull * svNull') + (dlFactor .* eye(nElem));

% % form array response matrix across fixed elevation
% UFixedEl = [sin(azAngles)*cos(zeros(K, 1)), cos(azAngles)*cos(zeros(K, 1)), sin(zeros(K, 1))].';
% S = exp(-1j * kc * (array.ElementPosition * UFixedEl));

% % solve for MVDR filter coefficients
% gain = max(S' * svSteer);
% wMVDR = (gain * R^-1 * svSteer) / (svSteer' * R^-1 * svSteer);

% % compute 2D array factor across azimuth sweep
% arrayFactorAz = S' * wMVDR;
% arrayFactorPerturbedAz = S' * (wMVDR .* svPhaseNoise .* svGainNoise);

% % form array response matrix across fixed azimuth
% UFixedAz = [sin(zeros(K, 1))*cos(elAngles), cos(zeros(K, 1))*cos(elAngles), sin(elAngles)].';
% S = exp(-1j * kc * (array.ElementPosition * UFixedAz));

% % compute 2D array factor across elevation sweep
% arrayFactorEl = S' * wMVDR;
% arrayFactorPerturbedEl = S' * (wMVDR .* svPhaseNoise .* svGainNoise);

% % plot 2d and 3d plots here



% %% Adaptive Beamforming with Coupling and Perturbations


% % define steering angle and form steering vector
% azSteer = azAngles(0);
% u = [sin(azSteer)*cos(0), cos(azSteer)*cos(0), sin(0)].';
% svSteer = exp(-1j * kc * (array.ElementPosition * u));

% % define null angle and form nulling vector
% azNull = azAngles(35);
% u = [sin(azNull)*cos(0), cos(azNull)*cos(0), sin(0)].';
% svNull = exp(-1j * kc * (array.ElementPosition * u));

% % form interference covariance matrix (and diagonally load)
% dlFactor = 10^(-30/20);
% R = (svNull * svNull') + (dlFactor .* eye(nElem));

% % form array response matrix across fixed elevation
% UFixedEl = [sin(azAngles)*cos(zeros(K, 1)), cos(azAngles)*cos(zeros(K, 1)), sin(zeros(K, 1))].';
% S = exp(-1j * kc * (array.ElementPosition * UFixedEl));

% % solve for MVDR filter coefficients
% gain = max(S' * svSteer);
% wMVDR = (gain * R^-1 * svSteer) / (svSteer' * R^-1 * svSteer);

% % compute 2D array factor across azimuth sweep
% arrayFactorAz = S' * wMVDR;
% arrayFactorCoupledAz = S' * (wMVDR .* ((eye(nElem) + sObj.Parameters) * ones(nElem, 1)) .* svPhaseNoise .* svGainNoise);

% % form array response matrix across fixed azimuth
% UFixedAz = [sin(zeros(K, 1))*cos(elAngles), cos(zeros(K, 1))*cos(elAngles), sin(elAngles)].';
% S = exp(-1j * kc * (array.ElementPosition * UFixedAz));

% % compute 2D array factor across elevation sweep
% arrayFactorEl = S' * wMVDR;
% arrayFactorCoupledEl = S' * (wMVDR .* ((eye(nElem) + sObj.Parameters) * ones(nElem, 1)) .* svPhaseNoise .* svGainNoise);

% % plot 2D and 3D patterns



% %% Plot patterns with coupling and perturbations