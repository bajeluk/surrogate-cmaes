% Array
% Submerged spheres, 3-tether arrangement, individual PTOs for each tether
% - fixed water depth
% - fixed submergence depth
% - one wave angle
% - fixed array layout
% - fixed optimal PTO parameters, detected for Tp = 9 sec, Hs = 2 m

% 26/10/2015    q-factor is changed to the RCW

function [RCWWeighted, ParrayAverage] = arrayBuoyPlacement_v2(radiusSphere, xSphere, ySphere, numFreq)

% Check whether we have enough inputs
% if length(radiusSphere) ~= numRows*numBodyInRow
%     error('number of specified radii is not equal to number of buoys in array')
% end

% Load file with optimal PTO parameters, maximum power absorption for the
% spheres of 5 different radii 5 m, 4 m, 3.2 m, 2.5 m, 2 m
load('LookUpTableForArray');

% Load file with Sydney probability data to calculate weighted q-factor and
% power absorption
load('rawDataSydney.mat');

% Water parameters
ro = 1020;              % kg/m^3, water density
g = 9.81;               % m/s^2
beta = 0;               % rad, wave propagation angle
subDepth = 6;           % m, submergence depth
oceanDepth = 50;        % m, ocean depth

numApprox = 4;      % number of equations

% Range of frequencies used to calculate array performance
[counts, w_array] = hist(fp, numFreq);

% Correct probability to match integral over freqencies to 1
countsTotal = sum(counts);
freqProbab = counts/countsTotal;


%% Array layout (layout #1)
% Calculation of the tether length depending on the water depth
alpha = acos(sqrt(1/3));            % tether inclination angle, instead of 54.7 deg
% tetherZ = (oceanDepth - subDepth);  % z-projection of the tether, from the sphere centre
% tetherL = tetherZ*sqrt(3);          % tether length, from the sphere centre (instead of (h-f)/cos(alpha))
% tetherXY = tetherZ*sqrt(2);         % x-y-projection of the sphere, instead of tetherL*sin(alpha) - it is also a radius of the circle for the hexagon
% 
% dY = tetherXY*sqrt(3);
% 
% xArray = zeros(numBodyInRow, numRows);
% yArray = zeros(numBodyInRow, numRows);
% 
% for row = 1:numRows
%     xArray(:, row) = sqrt(3)/2*(row-1)*ones(1,numBodyInRow);
%     if mod(row, 2) == 1
%         yArray(:, row) = 0:1:numBodyInRow-1;
%     else
%         yArray(:, row) = 1/2+(0:1:numBodyInRow-1);
%     end
% end
% 
numSphere = length(xSphere);
% xSphere = reshape(xArray, [numSphere, 1])*dY;
% ySphere = reshape(yArray, [numSphere, 1])*dY;
zSphere = -subDepth*ones(numSphere, 1);

%%
% Array structure
array.number = numSphere;
array.radius = radiusSphere;
array.sphereCoordinate(1,:) = xSphere;
array.sphereCoordinate(2,:) = ySphere;
array.sphereCoordinate(3,:) = zSphere;

% Wave structure
wave.waterDensity = ro;
wave.angle = beta;

% PTO specification
I3 = eye(3);

% Identification of vectors according to Scruggs (2013)
% Unit vectors along tethers (from ocean bottom to the attachment point)
es01 = [-sin(alpha);	0;                      cos(alpha)];
es02 = [sin(alpha)/2;	sqrt(3)/2*sin(alpha);	cos(alpha)];
es03 = [sin(alpha)/2;	-sqrt(3)/2*sin(alpha);	cos(alpha)];

% Kt and Ct from Scruggs (2013)
% Eq.11
Gt1 = (-I3*es01);
Gt2 = (-I3*es02);
Gt3 = (-I3*es03);

% Parameters for each sphere according to its radius
CtArray = zeros(numSphere*3);
KtArray = zeros(numSphere*3);
M = zeros(numSphere*3);         % Array mass matrix
Pi = zeros(numSphere, length(w_array));       % Absorbed power of individual buoys

for ii = 1:numSphere

    idx = 3*(ii-1)+(1:1:3);     % for 3 modes

    % Settings for the sphere of this radius
    a = radiusSphere(ii);

    fname = ['radius', num2str(a*10)];

    V = lookUpTable.(fname).volume;
    mass = lookUpTable.(fname).mass;
    kPTO = lookUpTable.(fname).kPTO;
    dPTO = lookUpTable.(fname).dPTO;

    % Net tether force
    Ft = (ro*g*V-mass*g);

    % Initial tether length
    s0n = (oceanDepth - subDepth - a*cos(alpha))/cos(alpha);

    % Initial tention in a leg
    t0 = Ft/(3*cos(alpha));
    gam0 = t0/s0n;

    % Eq.13
    Ct1 = dPTO*(Gt1*Gt1');
    Ct2 = dPTO*(Gt2*Gt2');
    Ct3 = dPTO*(Gt3*Gt3');

    % Eq.14
    Kt1 = ((kPTO-gam0)*(Gt1*Gt1') + gam0*I3);
    Kt2 = ((kPTO-gam0)*(Gt2*Gt2') + gam0*I3);
    Kt3 = ((kPTO-gam0)*(Gt3*Gt3') + gam0*I3);

    Ct = Ct1 + Ct2 + Ct3;
    Kt = Kt1 + Kt2 + Kt3;

    % PTO impedance for the array
    CtArray(idx, idx) = Ct;
    KtArray(idx, idx) = Kt;

    % Mass matrix
    M(idx, idx) = mass*I3;

    % Find power absorption of this sphere
    Pi100 = lookUpTable.(fname).power;
    
    % Interpolate power according to the new set of frequencies
    if length(w_array) ~= 100
        PiInterp = interp1(lookUpTable.radius50.freq, Pi100, w_array);
    else
        PiInterp = Pi100;
    end
    
    Pi(ii, :) = PiInterp;
end

Parray = zeros(1, length(w_array));
RCW = zeros(1, length(w_array));

% parpool();

parfor count_w = 1:length(w_array)
    
    w = w_array(count_w);
    K = w^2/g;

    % Wave power (infinite water depth)
    Pw = ro*g^2/(4*w);
    
    % Array hydrodynamics
    [A, B, X] = arraySubmergedSphereParfor(array, wave, w, K, numApprox, 1);
    
    % PTO impedance
    Zpto = (CtArray - 1i*KtArray/w);

    % Buoy impedance
    Zbuoy = 1i*(M + A)*w + B;

    % Velocity vector
    Uarray = (Zbuoy + Zpto)^(-1)*X;

    % Total absorbed power of the array
    Parray(count_w) = real(1/4*(Uarray'*X + X'*Uarray) - 1/2*Uarray'*B*Uarray);
    
    % RCW
    RCW(count_w) = Parray(count_w)/(2*sum(radiusSphere)*Pw);   % P/2a
end

% save('results', Parray, RCW);

% Weighted RCW according to the sea state probability
RCWWeighted = RCW*transpose(freqProbab)/sum(freqProbab);

% Average absorbed power of the array for this sea state probability
ParrayAverage = Parray*transpose(freqProbab)/sum(freqProbab);

end