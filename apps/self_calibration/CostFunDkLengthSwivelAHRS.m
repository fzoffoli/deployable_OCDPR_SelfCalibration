function [F,J] = CostFunDkLengthSwivelAHRS(cdpr_v, cdpr_p, measures, pose)
% Computes the residual vector and the jacobian for a
% state estimation using cable lengths, swivel angles and ahrs meas.

% sensor errors for normalization
swivel_noise = deg2rad(0.5);    %[rad]
length_noise = 0.02;             %[m]
AHRS_noise = deg2rad(0.01);         %[rad]


cable_lengths_meas  = measures.lengths;
swivel_angles_meas  = measures.swivels;
epsilon_meas        = measures.epsilon;

cdpr_v = UpdateIKZeroOrd(pose(1:3),pose(4:6),cdpr_p,cdpr_v);

cable_lengths = zeros(cdpr_p.n_cables,1);
swivel_angles = zeros(cdpr_p.n_cables,1);
for i = 1:cdpr_p.n_cables
    cable_lengths(i) = cdpr_v.cable(i).complete_length;
    swivel_angles(i) = cdpr_v.cable(i).swivel_ang;
end

F = diag([ones(cdpr_p.n_cables,1)./length_noise;ones(cdpr_p.n_cables,1)./swivel_noise; ...
            ones(3,1)./AHRS_noise])*[cable_lengths - cable_lengths_meas;
    swivel_angles - swivel_angles_meas;
    pose(4:6) - epsilon_meas];
J = diag([ones(cdpr_p.n_cables,1)./length_noise;ones(cdpr_p.n_cables,1)./swivel_noise; ...
            ones(3,1)./AHRS_noise])*[cdpr_v.analitic_jacobian_l';
    cdpr_v.analitic_jacobian_s';
    zeros(3), eye(3)];
end