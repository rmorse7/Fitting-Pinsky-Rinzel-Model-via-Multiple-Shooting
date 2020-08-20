function intervals = get_intervals(tspan, v)

starts = [];
ends = [];
intervals = [];

threshold = 0;
tn = length(tspan);
for i = 2:(tn-1)
    if v(i) < threshold && v(i+1) > threshold
        %starts = [starts, tspan(i)];
        intervals = [intervals, tspan(i)];%tspan(i)];
    elseif v(i) < v(i-1) && v(i) < v(i+1)
        %starts = [starts, tspan(i)];
        %intervals = [intervals, tspan(i)];%tspan(i)];
    elseif v(i) > threshold && v(i+1) < threshold
        %ends = [ends, tspan(i+1)];
        %intervals = [intervals, tspan(i+1)];%tspan(i+1)];
    elseif v(i) > v(i-1) && v(i) > v(i+1)
        %ends = [ends, tspan(i)];
        %intervals = [intervals, tspan(i)];
    end
end

%intervals = [starts;ends];
end