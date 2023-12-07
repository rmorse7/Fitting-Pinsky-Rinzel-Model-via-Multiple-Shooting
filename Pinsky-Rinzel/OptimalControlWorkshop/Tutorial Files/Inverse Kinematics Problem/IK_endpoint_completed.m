%=========================================================================%
% Endpoint function for the inverse kinematics problem
%=========================================================================%

function output = IK_endpoint_completed(input)

output.objective = input.phase.integral;

end