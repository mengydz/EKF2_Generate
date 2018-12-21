%% Set initial conditions
clear all;
dt = 1/100;
duration = 200;
indexLimit = round(duration/dt);
t = 0:dt:(duration-dt);
statesLog = zeros(9,indexLimit);
quatLog   = zeros(4,indexLimit);
velInnovLog = zeros(3,indexLimit);
decInnovLog = zeros(1,indexLimit);
velInnovVarLog = velInnovLog;
decInnovVarLog = decInnovLog;

eulerLog = zeros(3,indexLimit);
% Use a random initial orientation
yaw321 = 0.5;
eulerTruth = [0.3 0.4 0.5];
% quatTruth = [rand;randn;randn;randn];
quatTruth = EulToQuat(eulerTruth);
quatLength = sqrt(quatTruth(1)^2 + quatTruth(2)^2 + quatTruth(3)^2 + quatTruth(4)^2);
quatTruth = quatTruth / quatLength;
TbnTruth = Quat2Tbn(quatTruth);
% initialise the filter to level
  quat = [1;0;0;0];
% quat = [rand;randn;randn;randn];
% quatlen = sqrt(quat(1)^2 + quat(2)^2 + quat(3)^2 + quat(4)^2);
% quat = quat / quatlen;
states = zeros(9,1);
Tbn = Quat2Tbn(quat);
% define the earths truth magnetic field
magEarthTruth = [0.3;0.1;-0.5];
% define the assumed declination using th etruth field plus location
% variation
measDec = atan2(magEarthTruth(2),magEarthTruth(1)) + 2*pi/180*randn;
% define the magnetometer bias errors
magMeasBias = 0.02*[randn;randn;randn];
% define the state covariances with the exception of the quaternion covariances
Sigma_velNED = 0.5; % 1 sigma uncertainty in horizontal velocity components
Sigma_dAngBias  = 1*pi/180*dt; % 1 Sigma uncertainty in delta angle bias
Sigma_angErr = 1; % 1 Sigma uncertainty in angular misalignment (rad)
covariance   = single(diag([Sigma_angErr*[1;1;1];Sigma_velNED*[1;1;1];Sigma_dAngBias*[1;1;1]].^2));
%% external movement disturbance

acc_dis = zeros(3,indexLimit);

for i = 1:indexLimit
   if  i > indexLimit /3 && i < 2*indexLimit/3
       acc_dis(1,i) = -5; 
   end
end

%% Main Loop
headingAligned=0;
measVel = [0;0;0];
for index = 1:indexLimit
    % synthesise IMU measurements
    angRate = 0.0*[randn;randn;randn] + [0.2;0.001;0.5];
    accel = 0.0*[randn;randn;randn] + transpose(TbnTruth)*([0;0;-9.81]+acc_dis(:,index));
    % predict states
    [quat, states, Tbn, delAng, delVel]  = PredictStates(quat,states,angRate,accel,dt);
    statesLog(:,index) = states;
    quatLog(:,index) = quat;
    % predict covariance matrix
    covariance  = PredictCovariance(delAng,delVel,quat,states,covariance,dt);
    % synthesise velocity measurements
    measVel = measVel + acc_dis(:,index)*dt + 0.1*[randn;randn;randn];
%     measVel = [0;0;0];
    % fuse velocity measurements
    [quat,states,angErr,covariance,velInnov,velInnovVar] = FuseVelocity(quat,states,covariance,measVel);
    velInnovLog(:,index) = velInnov;
    velInnovVarLog(:,index) = velInnovVar;
%     % synthesise magnetometer measurements adding sensor bias
%     magBody = transpose(TbnTruth)*magEarthTruth + magMeasBias;
%     % fuse magnetometer measurements
%     if (index > 20 && headingAligned==0 && angErr < 1e-4)
%         quat = AlignHeading(quat,magBody,measDec);
%         headingAligned = 1;
%     end
%     if (headingAligned == 1)
%         [quat,states,covariance,decInnov,decInnovVar] = FuseMagnetometer(quat,states,covariance,magBody,measDec,Tbn);
%         decInnovLog(:,index) = decInnov;
%         decInnovVarLog(:,index) = decInnovVar;
%     end
    [quat, states,covariance, decInnov,decInnovVar] = fuseYaw321( quat,states,covariance, yaw321, Tbn);
    decInnovLog(:,index) = decInnov;
    decInnovVarLog(:,index) = decInnovVar;
    eulerLog(:,index) =  QuatToEul(quat);
end
figure;
plot(t,eulerLog(1,:),'r',t,eulerLog(2,:),'g',t,eulerLog(3,:),'b');grid on;