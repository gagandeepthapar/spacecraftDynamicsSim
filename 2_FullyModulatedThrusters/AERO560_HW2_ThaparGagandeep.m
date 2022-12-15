%% Gagandeep Thapar; AERO 560 HW2

% housekeeping
clc;
clearvars;
close all;

%% Problem 2

%% givens
orbit.h = 129654.3;    % [kg/m2] angular momentum
orbit.ecc = 0;  % [~] eccentricity
orbit.raan = 0;    % [deg] raan
orbit.inc = 0;     % [deg] inc
orbit.omega = 0;   % [deg] arg of perigee
orbit.theta = 0;    % [deg] true anomaly
orbit.mu = 398600;

sat.mass = 1200;    % [kg]  mass of sat
sat.r = 1.5;        % [m] radius of sat
sat.h = 10;         % [m] height of sat

sat.J = sat.mass/12 * [sat.h^2 + 3*sat.r^2, 0, 0;
                       0, sat.h^2 + 3*sat.r^2, 0;
                       0, 0, 6*sat.r^2];    % principal inertial matrix

sat.des_e_lvlh = [0;0;0];
sat.des_n_lvlh = 1;
sat.des_q_lvlh = [sat.des_e_lvlh;sat.des_n_lvlh];
sat.zeta = 0.65;    % [~] Dampening Coefficient
sat.ts = 30;    % [sec] settling time



%% initial conditions
[orbit.R0, orbit.V0] = COES2STATE(orbit.h, orbit.ecc, orbit.inc, orbit.raan, orbit.omega, orbit.theta, 398600);
orbit.R = orbit.R0;
orbit.V = orbit.V0;
orbit.T = 2*pi*norm(orbit.R0)^(1.5)/sqrt(orbit.mu);
orbit.n = 2*pi/orbit.T;

orbit.state0 = [orbit.R;orbit.V];

sat.w0_lvlh = [0.5;-7.27;3.0]*10^(-5); % initial ang vel
sat.e0_lvlh = [0.5;0.5;0.5];   % initial quat
sat.n0_lvlh = 0.5;
sat.q0_lvlh = [sat.e0_lvlh;sat.n0_lvlh];
sat.euls0_lvlh = quat_eul([sat.e0_lvlh;sat.n0_lvlh]); % initial euler

sat.state0_lvlh = [sat.w0_lvlh;sat.e0_lvlh;sat.n0_lvlh;sat.euls0_lvlh];  % initial state

sat.wn = log(0.02*sqrt(1 - sat.zeta^2))/(-1*sat.zeta*sat.ts);
sat.zeta = sat.zeta*eye(3);
sat.wn = sat.wn*eye(3);

sat.Kd = 2*sat.J*sat.zeta*sat.wn;
sat.Kp = 2*(sat.J)*(sat.wn^2);

%% run sim
t_max = 3600;    % [sec] sim time
out_A = sim('hw2n_lawA', t_max);
out_B = sim('hw2n_lawB', t_max);
out_C = sim('hw2n_lawC', t_max);

%% unpack data: A
A_t = squeeze(out_A.tout)';

A_orbit.w_lvlh = squeeze(out_A.w_LVLH_ECI)';
A_orbit.q_lvlh = squeeze(out_A.q_LVLH_ECI)';
A_orbit.eul_lvlh = squeeze(out_A.eul_LVLH_ECI)';

A_orbit.dist_torque_lvlh = squeeze(out_A.dist_torque)';

A_sat.w_lvlh = squeeze(out_A.w_Body_LVLH)';
A_sat.q_lvlh = squeeze(out_A.q_Body_LVLH)';
A_sat.eul_lvlh = squeeze(out_A.eul_Body_LVLH)';

A_sat.w_eci = squeeze(out_A.w_Body_ECI)';
A_sat.q_eci = squeeze(out_A.q_Body_ECI)';
A_sat.eul_eci = squeeze(out_A.eul_Body_ECI)';

A_command_torque = squeeze(out_A.command_torque)';

%% unpack data: B
B_t = squeeze(out_B.tout)';

B_orbit.w_lvlh = squeeze(out_B.w_LVLH_ECI)';
B_orbit.q_lvlh = squeeze(out_B.q_LVLH_ECI)';
B_orbit.eul_lvlh = squeeze(out_B.eul_LVLH_ECI)';

B_orbit.dist_torque_lvlh = squeeze(out_B.dist_torque)';

B_sat.w_lvlh = squeeze(out_B.w_Body_LVLH)';
B_sat.q_lvlh = squeeze(out_B.q_Body_LVLH)';
B_sat.eul_lvlh = squeeze(out_B.eul_Body_LVLH)';

B_sat.w_eci = squeeze(out_B.w_Body_ECI)';
B_sat.q_eci = squeeze(out_B.q_Body_ECI)';
B_sat.eul_eci = squeeze(out_B.eul_Body_ECI)';

B_command_torque = squeeze(out_B.command_torque)';

%% unpack data: C
C_t = squeeze(out_C.tout)';

C_orbit.w_lvlh = squeeze(out_C.w_LVLH_ECI)';
C_orbit.q_lvlh = squeeze(out_C.q_LVLH_ECI)';
C_orbit.eul_lvlh = squeeze(out_C.eul_LVLH_ECI)';

C_orbit.dist_torque_lvlh = squeeze(out_C.dist_torque)';

C_sat.w_lvlh = squeeze(out_C.w_Body_LVLH)';
C_sat.q_lvlh = squeeze(out_C.q_Body_LVLH)';
C_sat.eul_lvlh = squeeze(out_C.eul_Body_LVLH)';

C_sat.w_eci = squeeze(out_C.w_Body_ECI)';
C_sat.q_eci = squeeze(out_C.q_Body_ECI)';
C_sat.eul_eci = squeeze(out_C.eul_Body_ECI)';

C_command_torque = squeeze(out_C.command_torque)';

%% plot data: A
figure

subplot(3,3,1)
hold on
plot(A_t, A_sat.q_lvlh(:,1))
plot(A_t, A_sat.q_lvlh(:,2))
plot(A_t, A_sat.q_lvlh(:,3))
plot(A_t, A_sat.q_lvlh(:,4))
hold off
grid on
title('Quaternion from Body to LVLH')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,2)
hold on
plot(A_t, A_orbit.q_lvlh(:,1))
plot(A_t, A_orbit.q_lvlh(:,2))
plot(A_t, A_orbit.q_lvlh(:,3))
plot(A_t, A_orbit.q_lvlh(:,4))
hold off
grid on
title('Quaternion from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,3)
hold on
plot(A_t, A_sat.q_eci(:,1))
plot(A_t, A_sat.q_eci(:,2))
plot(A_t, A_sat.q_eci(:,3))
plot(A_t, A_sat.q_eci(:,4))
hold off
grid on
title('Quaternion from Body to ECI')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,4)
hold on
plot(A_t, A_sat.eul_lvlh(:,1))
plot(A_t, A_sat.eul_lvlh(:,2))
plot(A_t, A_sat.eul_lvlh(:,3))
hold off
grid on
title('Euler Angles from Body to LVLH')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')

subplot(3,3,5)
hold on
plot(A_t, A_orbit.eul_lvlh(:,1))
plot(A_t, A_orbit.eul_lvlh(:,2))
plot(A_t, A_orbit.eul_lvlh(:,3))
hold off
grid on
title('Euler Angles from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')

subplot(3,3,6)
hold on
plot(A_t, A_sat.eul_eci(:,1))
plot(A_t, A_sat.eul_eci(:,2))
plot(A_t, A_sat.eul_eci(:,3))
hold off
grid on
title('Euler Angles from Body to ECI')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')


subplot(3,3,7)
hold on
plot(A_t, A_sat.w_lvlh(:,1))
plot(A_t, A_sat.w_lvlh(:,2))
plot(A_t, A_sat.w_lvlh(:,3))
hold off
grid on
title('\omega from Body to LVLH')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

subplot(3,3,8)
hold on
plot(A_t, A_orbit.w_lvlh(:,1))
plot(A_t, A_orbit.w_lvlh(:,2))
plot(A_t, A_orbit.w_lvlh(:,3))
hold off
grid on
title('\omega from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

subplot(3,3,9)
hold on
plot(A_t, A_sat.w_eci(:,1))
plot(A_t, A_sat.w_eci(:,2))
plot(A_t, A_sat.w_eci(:,3))
hold off
grid on
title('\omega from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

sgtitle('T_{C} = -K_psign(n_e)\epsilon_e - K_d\omega')

%% plot data: B
figure

subplot(3,3,1)
hold on
plot(B_t, B_sat.q_lvlh(:,1))
plot(B_t, B_sat.q_lvlh(:,2))
plot(B_t, B_sat.q_lvlh(:,3))
plot(B_t, B_sat.q_lvlh(:,4))
hold off
grid on
title('Quaternion from Body to LVLH')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,2)
hold on
plot(B_t, B_orbit.q_lvlh(:,1))
plot(B_t, B_orbit.q_lvlh(:,2))
plot(B_t, B_orbit.q_lvlh(:,3))
plot(B_t, B_orbit.q_lvlh(:,4))
hold off
grid on
title('Quaternion from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,3)
hold on
plot(B_t, B_sat.q_eci(:,1))
plot(B_t, B_sat.q_eci(:,2))
plot(B_t, B_sat.q_eci(:,3))
plot(B_t, B_sat.q_eci(:,4))
hold off
grid on
title('Quaternion from Body to ECI')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,4)
hold on
plot(B_t, B_sat.eul_lvlh(:,1))
plot(B_t, B_sat.eul_lvlh(:,2))
plot(B_t, B_sat.eul_lvlh(:,3))
hold off
grid on
title('Euler Angles from Body to LVLH')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')

subplot(3,3,5)
hold on
plot(B_t, B_orbit.eul_lvlh(:,1))
plot(B_t, B_orbit.eul_lvlh(:,2))
plot(B_t, B_orbit.eul_lvlh(:,3))
hold off
grid on
title('Euler Angles from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')

subplot(3,3,6)
hold on
plot(B_t, B_sat.eul_eci(:,1))
plot(B_t, B_sat.eul_eci(:,2))
plot(B_t, B_sat.eul_eci(:,3))
hold off
grid on
title('Euler Angles from Body to ECI')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')


subplot(3,3,7)
hold on
plot(B_t, B_sat.w_lvlh(:,1))
plot(B_t, B_sat.w_lvlh(:,2))
plot(B_t, B_sat.w_lvlh(:,3))
hold off
grid on
title('\omega from Body to LVLH')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

subplot(3,3,8)
hold on
plot(B_t, B_orbit.w_lvlh(:,1))
plot(B_t, B_orbit.w_lvlh(:,2))
plot(B_t, B_orbit.w_lvlh(:,3))
hold off
grid on
title('\omega from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

subplot(3,3,9)
hold on
plot(B_t, B_sat.w_eci(:,1))
plot(B_t, B_sat.w_eci(:,2))
plot(B_t, B_sat.w_eci(:,3))
hold off
grid on
title('\omega from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

sgtitle('T_{C} = -K_psign(n_e)\epsilon_e - K_d(1-\epsilon_e^T\epsilon_e)\omega')

%% plot data: C
figure

subplot(3,3,1)
hold on
plot(C_t, C_sat.q_lvlh(:,1))
plot(C_t, C_sat.q_lvlh(:,2))
plot(C_t, C_sat.q_lvlh(:,3))
plot(C_t, C_sat.q_lvlh(:,4))
hold off
grid on
title('Quaternion from Body to LVLH')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,2)
hold on
plot(C_t, C_orbit.q_lvlh(:,1))
plot(C_t, C_orbit.q_lvlh(:,2))
plot(C_t, C_orbit.q_lvlh(:,3))
plot(C_t, C_orbit.q_lvlh(:,4))
hold off
grid on
title('Quaternion from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,3)
hold on
plot(C_t, C_sat.q_eci(:,1))
plot(C_t, C_sat.q_eci(:,2))
plot(C_t, C_sat.q_eci(:,3))
plot(C_t, C_sat.q_eci(:,4))
hold off
grid on
title('Quaternion from Body to ECI')
xlabel('Time [sec]')
ylabel('Quaternion')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')

subplot(3,3,4)
hold on
plot(C_t, C_sat.eul_lvlh(:,1))
plot(C_t, C_sat.eul_lvlh(:,2))
plot(C_t, C_sat.eul_lvlh(:,3))
hold off
grid on
title('Euler Angles from Body to LVLH')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')

subplot(3,3,5)
hold on
plot(C_t, C_orbit.eul_lvlh(:,1))
plot(C_t, C_orbit.eul_lvlh(:,2))
plot(C_t, C_orbit.eul_lvlh(:,3))
hold off
grid on
title('Euler Angles from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')

subplot(3,3,6)
hold on
plot(C_t, C_sat.eul_eci(:,1))
plot(C_t, C_sat.eul_eci(:,2))
plot(C_t, C_sat.eul_eci(:,3))
hold off
grid on
title('Euler Angles from Body to ECI')
xlabel('Time [sec]')
ylabel('Euler Angle [rad]')
legend('\phi', '\theta', '\psi')


subplot(3,3,7)
hold on
plot(C_t, C_sat.w_lvlh(:,1))
plot(C_t, C_sat.w_lvlh(:,2))
plot(C_t, C_sat.w_lvlh(:,3))
hold off
grid on
title('\omega from Body to LVLH')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

subplot(3,3,8)
hold on
plot(C_t, C_orbit.w_lvlh(:,1))
plot(C_t, C_orbit.w_lvlh(:,2))
plot(C_t, C_orbit.w_lvlh(:,3))
hold off
grid on
title('\omega from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

subplot(3,3,9)
hold on
plot(C_t, C_sat.w_eci(:,1))
plot(C_t, C_sat.w_eci(:,2))
plot(C_t, C_sat.w_eci(:,3))
hold off
grid on
title('\omega from LVLH to ECI')
xlabel('Time [sec]')
ylabel('Rate [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z')

sgtitle('T_{C} = -K_psign(n_e)\epsilon_e - K_d(1+\epsilon_e^T\epsilon_e)\omega')

%% plot misc

for i = 1:length(A_command_torque)
    A_tq(i) = norm(A_command_torque(i,:));
end

for i = 1:length(B_command_torque)
    B_tq(i) = norm(B_command_torque(i,:));
end

for i = 1:length(C_command_torque)
    C_tq(i) = norm(C_command_torque(i,:));
end

figure
hold on
plot(A_t, A_tq)
plot(B_t, B_tq)
plot(C_t, C_tq)
hold off
title('Torque Exerted')
xlabel('Time [sec]')
ylabel('Norm of Torque [Nm]')
legend('Control Law A', 'Control Law B', 'Control Law C')

fprintf('~~~~~~~~~~~~~~~~~~~~\n')
fprintf('4. The satellite can maintain pointing with the disturbance torque. Setting the disturbance torque to 0 results in a similar result with the satellite reaching steady state in less time.\n')

fprintf('~~~~~~~~~~~~~~~~~~~~\n')
fprintf('5. The max torque exerted (in the second control law) is %.2f Nm. With a height of %.2f m, the satellite must exert %.2f N of thrust. No EP is currently capable of producing that much thrust.\n', max(B_tq), sat.h, max(B_tq)/(sat.h/2));

fprintf('~~~~~~~~~~~~~~~~~~~~\n')
fprintf('6. Each control law requires different amount of torque for different amounts of time, however, each reach steady state in a similar amount of time. Control Law B has a greater peak than its counterparts as expected due to the control law looking at the difference in the quaternion near equilibrium (which will cause a larger torque requirement initially).\n');

fprintf('~~~~~~~~~~~~~~~~~~~~\n')

function eul = quat_eul(q)

    n = q(4);
    e = q(1:3);

    q = [n, e(1), e(2), e(3)];

    phi = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2));
    theta = asin(2*(q(1)*q(3) - q(4)*q(2)));
    psi = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1-2*(q(3)^2 + q(4)^2));

    eul = [phi; theta; psi];

end

