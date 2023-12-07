function output = squat_endpoint(input)

% Extrack initial and final time and states
state1_final = input.phase(1).finalstate;
state2_initial = input.phase(2).initialstate;
t1_final = input.phase(1).finaltime;
t2_initial = input.phase(2).initialtime;

% Output objective and event constraint values
q = input.phase(1).integral+input.phase(2).integral;
output.objective = q;

output.eventgroup.event = [state1_final - state2_initial,t1_final - t2_initial];
end