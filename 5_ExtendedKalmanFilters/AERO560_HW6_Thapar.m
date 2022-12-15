%% Gagandeep Thapar; HW6 AERO 560

%% Housekeeping

clc;
clear all;
close all;

% random seed
seed = floor(sum(clock));
%% Part 1

% givens
p1.wn = 2;
p1.zeta = 0.7;
p1.delT = 0.1;

p1.x0 = 0.05;
p1.xd0 = -0.01;
p1.Q = 0.1;
p1.R = 0.1;
p1.sigma = 0.0035;
p1.u = 0;

% work
p1.A = [0 1; -p1.wn^2 -2*p1.zeta*p1.wn];
p1.B = [0;p1.wn^2];
p1.C = eye(2);
p1.D = zeros(2,1);
p1.state0 = [p1.x0;p1.xd0];

p1.sysC = ss(p1.A, p1.B, p1.C, p1.D);
p1.sysD = c2d(p1.sysC, p1.delT);

p1.F = p1.sysD.A;
p1.G = p1.sysD.B;
p1.H = p1.sysD.C;

p1.L = p1.G;
p1.M = p1.H;

% simulate
p1.outC = sim("AERO560_HW6Sim_Thapar.slx");
p1.t = squeeze(p1.outC.tout)';
p1.yC = squeeze(p1.outC.y_Clean)';
p1.yN = squeeze(p1.outC.y_Noisy)';
p1.xHat = squeeze(p1.outC.xHat);

% plot
figure
subplot(2,1,1)
hold on
plot(p1.t, p1.yC(:,1))
plot(p1.t, p1.yN(:,1))
plot(p1.t, p1.xHat(:,1),'k', 'LineWidth', 2)
hold off
legend('Clean', 'Noisy', 'Kalman Estimate')
title('Position Estimate')
xlabel('Time [sec]')
ylabel('Position [m]')

subplot(2,1,2)
hold on
plot(p1.t, p1.yC(:,2))
plot(p1.t, p1.yN(:,2))
plot(p1.t, p1.xHat(:,2),'k', 'LineWidth', 2)
hold off
legend('Clean', 'Noisy', 'Kalman Estimate')
title('Velocity Estimate')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

%% Part 2

% givens
p2.delT = 1;
p2.sA = [1;0;0];
p2.bA = [0.3815;-0.0969;0.9193];

p2.u = 0;
p2.J = diag([27, 17, 25]);

p2.w0 = [0.05;-0.05;0.05];
p2.e0 = sin(0.5)/sqrt(3) * [1;1;1];
p2.n0 = cos(0.5);

p2.wHat0 = zeros(3,1);
p2.eHat0 = zeros(3,1);
p2.nHat0 = 1;

p2.Pk0 = 0.1*eye(7);

p2.sigQ = 0.001;
p2.Q = p2.sigQ^2;

p2.sigM = 0.01;
p2.sigS = 0.005;
p2.R = diag([p2.sigM^2, p2.sigM^2, p2.sigM^2, p2.sigS^2, p2.sigS^2, p2.sigS^2]);

p2.v = 0;

% simulate
p2.out = sim('AERO560_HW6_Sim2_Thapar.slx', 100);

p2.t = squeeze(p2.out.tout);
p2.wTrue = squeeze(p2.out.wTrue);
p2.eTrue = squeeze(p2.out.eTrue);
p2.nTrue = squeeze(p2.out.nTrue)';
p2.xHat = squeeze(p2.out.xHat);

p2.wHat = p2.xHat(:,1:3);
p2.eHat = p2.xHat(:,4:6);
p2.nHat = p2.xHat(:,7);

% plot
figure
subplot(4,1,1)
hold on
plot(p2.t, p2.eTrue(:,1), '--')
plot(p2.t, p2.eHat(:,1), '-o')
hold off
xlabel('Time [sec]')
ylabel('\epsilon_1')
title('\epsilon_1')
legend('True', 'Kalman Estimate')

subplot(4,1,2)
hold on
plot(p2.t, p2.eTrue(:,2),'--')
plot(p2.t, p2.eHat(:,2), '-o')
hold off
xlabel('Time [sec]')
ylabel('\epsilon_2')
title('\epsilon_2')
legend('True', 'Kalman Estimate')

subplot(4,1,3)
hold on
plot(p2.t, p2.eTrue(:,3),'--')
plot(p2.t, p2.eHat(:,3), '-o')
hold off
xlabel('Time [sec]')
ylabel('\epsilon_3')
title('\epsilon_3')
legend('True', 'Kalman Estimate')

subplot(4,1,4)
hold on
plot(p2.t, p2.nTrue, '--')
plot(p2.t, p2.nHat, '-o')
hold off
xlabel('Time [sec]')
ylabel('\eta')
title('\eta')
legend('True', 'Kalman Estimate')
sgtitle('Estimate of Attitude Over Time')

figure
subplot(3,1,1)
hold on
plot(p2.t, p2.wTrue(:,1), '--')
plot(p2.t, p2.wHat(:,1), '-o')
hold off
xlabel('Time [sec]')
ylabel('\omega_x [rad/s]')
title('\omega_x')
legend('True', 'Kalman Estimate')

subplot(3,1,2)
hold on
plot(p2.t, p2.wTrue(:,2), '--')
plot(p2.t, p2.wHat(:,2), '-o')
hold off
xlabel('Time [sec]')
ylabel('\omega_y [rad/s]')
title('\omega_y')
legend('True', 'Kalman Estimate')

subplot(3,1,3)
hold on
plot(p2.t, p2.wTrue(:,3), '--')
plot(p2.t, p2.wHat(:,3), '-o')
hold off
xlabel('Time [sec]')
ylabel('\omega_z [rad/s]')
title('\omega_z')
legend('True', 'Kalman Estimate')
sgtitle('Estimate of Angular Velocity Over Time')

figure
hold on
plot(p2.t, p2.eTrue(:,1), '--r')
plot(p2.t, p2.eTrue(:,2), '--g')
plot(p2.t, p2.eTrue(:,3), '--b')
plot(p2.t, p2.nTrue, '--m')

plot(p2.t, p2.eHat(:,1), '-or');
plot(p2.t, p2.eHat(:,2), '-og');
plot(p2.t, p2.eHat(:,3), '-ob');
plot(p2.t, p2.nHat, '-om');

hold off
legend('\epsilon_1', '\epsilon_2', '\epsilon_3', '\eta')
xlabel('Time [sec]')
ylabel('Quaternion [~]')
title('Attitude Estimate')

figure
hold on
plot(p2.t, p2.wTrue(:,1), '--r')
plot(p2.t, p2.wTrue(:,2), '--g')
plot(p2.t, p2.wTrue(:,3), '--b')

plot(p2.t, p2.wHat(:,1), '-or')
plot(p2.t, p2.wHat(:,2), '-og')
plot(p2.t, p2.wHat(:,3), '-ob')

hold off
legend('\omega_x', '\omega_y', '\omega_z')
xlabel('Time [sec]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity Estimate')

