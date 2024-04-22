mass = 68.1;          
g = 9.8;            
drag_coefficient = 0.25; 
time_step = 0.1;      
total_time = 12; 

velocity = 0;         
position = 0;         

terminal_velocity = sqrt((mass * g) / drag_coefficient);

time_values = 0:time_step:total_time;
velocity_values = zeros(size(time_values));
position_values = zeros(size(time_values));

% Euler method
for i = 1:length(time_values)

    drag_force = drag_coefficient * velocity^2;

    acceleration = (mass * g - drag_force) / mass;

    velocity = velocity + acceleration * time_step;
    position = position + velocity * time_step;

    velocity_values(i) = velocity;
    position_values(i) = position;
end

figure;
subplot(2,1,1);
plot(time_values, velocity_values, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs Time');

subplot(2,1,2);
plot(time_values, position_values, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position (m)');
title('Position vs Time');

disp(['Terminal Velocity: ' num2str(terminal_velocity) ' m/s']);