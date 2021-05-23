% plot( [-0.05, 0.05, 0.05, -0.05, -0.05] , [-0.05,-0.05,0.2,0.2,-0.05] , 'color' , [72/255 130/255 197/255], 'LineWidth', 2.0);
% hold on;
clear;
Ts = 0.02;  % Sample Time
N = 5;    % Prediction horizon
Duration = 1.0; % Simulation horizon
number = 100;

% For a usual control 0.1 , 5, 20

% Things specfic to the model starts here 

% Limits of initial [-1,1]^2

global simulation_result;
global disturb_range;

disturb_range = 0; % Disturbance range

formatSpec = '%f %f\n';

fileID = fopen('nn_1_sigmoid.txt','w');

% figure;
for m=1:number
% Setting the initial state
disp(num2str(m) + "-th simulation starts: ")

x0 = 0.09 + 0.02*rand(1);
x1 = - 0.21 + 0.02*rand(1);
x2 = - 0.21 + 0.02*rand(1);
x3 = 0.29 + 0.02*rand(1);
x4 = 0.04 + 0.02*rand(1);
x5 = - 0.01 + 0.02*rand(1);

x = [x0;x1;x2;x3;x4;x5];

simulation_result = x;

% Things specfic to the model  ends here

options = optimoptions('fmincon','Algorithm','sqp','Display','none');
uopt = zeros(N,1);

u_max = 0;

% Apply the control input constraints
LB = -2*ones(N,1);
UB = 2*ones(N,1);

x_now = zeros(2,1);
x_next = zeros(2,1);

x_now = x;

% Start simulation
trace = [x_now];
for ct = 1:(Duration/Ts)
    
     u = NN_output(x_now,'CLF_controller_layer_num_3');
    % disp(u);
      
%     x_next = system_eq_NN(x_now, Ts, u);
      x_next = system_eq_dis(x_now, Ts, u);

      x = x_next;
      x_now = x_next;
      trace = [trace x_now];
% Save plant states for display.
%    xHistory = [xHistory x]; 


%     Printing stuff

%    disp('x(1) = '); x(1)
%    disp('x(2) = '); x(2)

%    disp('Control input = '); u
    
        
end

% fclose(file);

% figure;
plot(simulation_result(5,:),simulation_result(6,:), 'color', [210/255, 95/255, 95/255]);
% title('Benchmark 1 (sigmoid)', 'FontSize', 14)
% xlabel('x1', 'FontSize', 14);
% ylabel('x2', 'FontSize', 14);
set(gca,'FontSize',16)
hold on;

% fprintf(fileID, formatSpec, trace(1:2,:));
% fprintf(fileID, '%c\n', ['*']);

end

% fclose(fileID);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'../Benchmarks/benchmark1_sigmoid','-dpdf')
% export_fig ../Benchmarks/benchmark1_sigmoid.pdf