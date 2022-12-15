%% Gagandeep Thapar
% AERO 560; HW 4

%% Housekeeping

clc;
clear all;
close all;


%% Givens and Setup

% problem 2
trigger.delOff = 0.1;
trigger.delOn = 0.2;
trigger.mag = 1;

trigger2.delOff = 0.001;
trigger2.delOn = 0.002;

sat.J = 100;

sat.p1.command = 0;
sat.p1.initial = 5;
sat.p1.tau = 5;



% problem 3
sat.mass = 750;
sat.side = 1;
sat.dx = sat.side/2;
sat.dy = sat.side/2;
sat.dz = sat.side/2;

sat.p3.J = sat.mass * sat.side^2 / 6 * eye(3);


sat.zeta = 0.65;    % [~] Dampening Coefficient
sat.ts = 30;    % [sec] settling time
sat.wn = log(0.02*sqrt(1 - sat.zeta^2))/(-1*sat.zeta*sat.ts);
sat.zeta = sat.zeta*eye(3);
sat.wn = sat.wn*eye(3);

sat.Kd = 2*sat.p3.J*sat.zeta*sat.wn;
sat.Kp = 2*(sat.p3.J)*(sat.wn^2);

sat.T = 50;
sat.HF = sat.T * [-1 0 0 -1 0 0 1 0 0 1 0 0;
          0 1 0 0 -1 0 0 1 0 0 -1 0;
          0 0 1 0 0 -1 0 0 -1 0 0 1];
sat.HM = [0 sat.dz -sat.dy 0 sat.dz -sat.dy 0 -sat.dz sat.dy 0 -sat.dz sat.dy;
          sat.dz 0 -sat.dx -sat.dz 0 sat.dx sat.dz 0 -sat.dx -sat.dz 0 sat.dx;
          -sat.dy sat.dx 0 sat.dy -sat.dx 0 sat.dy -sat.dx 0 -sat.dy sat.dx 0];

sat.H = [sat.HF;sat.HM];

sat.p3.psi0 = 23;
sat.p3.theta0 = -16;
sat.p3.phi0 = 42;

sat.p3.euls0_lvlh = [sat.p3.psi0;sat.p3.theta0;sat.p3.phi0];
sat.p3.w0_lvlh = [0;0;0];

sat.p3.quat0 = eul_quat(sat.p3.phi0, sat.p3.theta0, sat.p3.psi0);
sat.p3.des_quat = [0;0;0;1];

%% Simulation
p1 = sim('HW4Model.slx', 1000);
p3_sim = sim('HW4ThrusterModel.slx', 100);

%% Unpack data
% P2
time = p1.tout;

sat.switch.theta = p1.switch_theta;
sat.switch.thetaDot = p1.switch_thetaDot;

sat.schmitt.theta = p1.schmitt_theta;
sat.schmitt.thetaDot = p1.schmitt_thetaDot;

% P3
p3.time = squeeze(p3_sim.tout)';
p3.req_torque = squeeze(p3_sim.command_torque)';
p3.sat.w = squeeze(p3_sim.w_Body_LVLH)';
p3.sat.eul = squeeze(p3_sim.eul_Body_LVLH)';
p3.sat.q = squeeze(p3_sim.q_Body_LVLH)';

p3.allocation = squeeze(p3_sim.thrust_alloc)';
p3.thrust_sol = squeeze(p3_sim.thrust_sol)';

%% Plot Data

% switch switching lines
x = [-5 5];
tau_line = -1/sat.p1.tau * x;

figure
hold on
plot(sat.switch.theta, sat.switch.thetaDot);
plot(x, tau_line, 'r--')
scatter(0,0,'ro')
hold off
xlabel('Theta')
ylabel('ThetaDot')
title('Switch Control: Theta VS ThetaDot')
legend('State Space', 'Switching Line', 'Set Point')
axis padded


figure
hold on
plot(sat.schmitt.theta, sat.schmitt.thetaDot);
% plot(x, tau_line, 'r--')
scatter(0,0,'ro')
hold off
xlabel('Theta')
ylabel('ThetaDot')
title('Schmitt Trigger: Theta VS ThetaDot')
axis padded
% legend('State Space', 'Switching Line', 'Set Point')

figure
hold on
plot(time, sat.switch.theta, 'r')
plot(time, sat.schmitt.theta, 'b')
hold off
xlabel('Time [sec]')
ylabel('Theta')
title('Time VS Theta')
axis padded
legend('Switch Control', 'Schmitt Trigger')

%% Plot Data (P3)
% close all;

figure
subplot(3,1,1)
plot(p3.time, p3.sat.q);
xlabel('Time [s]')
ylabel('Quaternion')
title('Body/LVLH Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,1,2)
plot(p3.time, p3.sat.eul);
xlabel('Time [s]')
ylabel('Euler Angle')
title('Euler Angles [deg]')
legend('Phi', 'Theta', 'Psi')

subplot(3,1,3)
plot(p3.time, p3.sat.w);
xlabel('Time [s]')
ylabel('Angular Rates [deg/s]')
title('Angular Velocity')
legend('\omega_1', '\omega_2', '\omega_3')

figure
subplot(3,1,1)
plot(p3.time, p3.req_torque);
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('Required Torque')

subplot(3,1,2)
plot(p3.time, p3.thrust_sol);
xlabel('Time [s]')
ylabel('Thrust Allocation [Nm]')
title('Thrust Allocation')

subplot(3,1,3)
plot(p3.time, p3.req_torque-p3.thrust_sol);
xlabel('Time [s]')
ylabel('Error [Nm]')
title('Thrust Allocation/Requirement Difference')

figure
[~,a] = size(p3.allocation);
hold on
for i = 1:a
p3.allocation(:,i) = p3.allocation(:,i)/sat.T;
plot(p3.time, p3.allocation(:,i));
end
hold off
xlabel('Time [s]')
ylabel('Thruster Allocation [Per Thruster, Ratio]')
%% functions
function eul = quat_eul(q)

    n = q(4);
    e = q(1:3);

    q = [n, e(1), e(2), e(3)];

    phi = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2));
    theta = asin(2*(q(1)*q(3) - q(4)*q(2)));
    psi = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1-2*(q(3)^2 + q(4)^2));

    eul = [phi; theta; psi];

end

function quat = eul_quat(phi, theta, psi)

    function Q = Qx(theta)
        Q = [1 0 0;
            0 cosd(theta) sind(theta);
            0 -sind(theta) cosd(theta)];
    end

    function Q = Qy(theta)
        Q = [cosd(theta) 0 -sind(theta);
                0 1 0;
                sind(theta) 0 cosd(theta)];
    end
    
    function Q = Qz(theta)
        Q = [cosd(theta) sind(theta) 0;
                -sind(theta) cosd(theta) 0;
                0 0 1];
    end

    C = Qx(phi)*Qy(theta)*Qz(psi);

    e = zeros(3,1);

    eta = 0.5 * sqrt(1 + trace(C));
    e(1) = 0.25 * (C(2,3) - C(3,2))/eta;
    e(2) = 0.25 * (C(3,1) - C(1,3))/eta;
    e(3) = 0.25 * (C(1,2) - C(2,1))/eta;

    quat = [e;eta];

end