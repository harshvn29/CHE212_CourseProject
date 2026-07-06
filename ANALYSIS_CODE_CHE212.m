
% 2D HEAT TRANSFER : FINAL PRESENTATION MODEL
% Logic: Conduction + Convection (h=10) + Radiation + Face Loss


clear; clc; close all;

%% 1. PARAMETERS & CALIBRATED CONSTANTS
N = 40;                 
L = 0.20;               
dx = L / (N - 1);

Th = 63;                
Tinf = 33;              
Ts_avg = 46;            

k = 15;                 
t = 0.004;              
h_conv = 10;            
k_ins = 0.045;          
L_ins = 0.004;          
epsi = 0.90;            
sigma = 5.67e-8;        

Tsur_K = Tinf + 273.15; 
m_sq = (2 * k_ins) / (k * t * L_ins);
loss_f = m_sq * dx^2;

%% 2. INITIALIZATION
T = ones(N, N) * ((Th + Tinf)/2); 
T(N, :) = Th;           

%% 3. ITERATIVE SOLVER (Gauss-Seidel)
fprintf('Simulating steady-state with all loss mechanisms...\n');
for iter = 1:200000
    T_old = T;
    
    for i = 2:N-1
        for j = 2:N-1
            T(i,j) = (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) + loss_f*Ts_avg) / (4 + loss_f);
        end
    end
    
    for i = 2:N-1
        TkL = T(i, 1) + 273.15;
        h_radL = epsi * sigma * (TkL + Tsur_K) * (TkL^2 + Tsur_K^2);
        h_effL = h_conv + h_radL;
        T(i, 1) = (0.5*k*(T(i+1, 1) + T(i-1, 1)) + k*T(i, 2) + h_effL*dx*Tinf + 0.5*k*loss_f*Ts_avg) / (2*k + h_effL*dx + 0.5*k*loss_f);
              
        TkR = T(i, N) + 273.15;
        h_radR = epsi * sigma * (TkR + Tsur_K) * (TkR^2 + Tsur_K^2);
        h_effR = h_conv + h_radR;
        T(i, N) = (0.5*k*(T(i+1, N) + T(i-1, N)) + k*T(i, N-1) + h_effR*dx*Tinf + 0.5*k*loss_f*Ts_avg) / (2*k + h_effR*dx + 0.5*k*loss_f);
    end

    for j = 2:N-1
        TkB = T(1, j) + 273.15;
        h_radB = epsi * sigma * (TkB + Tsur_K) * (TkB^2 + Tsur_K^2);
        h_effB = h_conv + h_radB;
        T(1, j) = (0.5*k*(T(1, j+1) + T(1, j-1)) + k*T(2, j) + h_effB*dx*Tinf + 0.5*k*loss_f*Ts_avg) / (2*k + h_effB*dx + 0.5*k*loss_f);
    end
    
    T(1, 1) = (k*(T(2, 1) + T(1, 2)) + 2*(h_conv+5)*dx*Tinf + 0.25*k*loss_f*Ts_avg) / (2*k + 2*(h_conv+5)*dx + 0.25*k*loss_f);
    T(1, N) = (k*(T(2, N) + T(1, N-1)) + 2*(h_conv+5)*dx*Tinf + 0.25*k*loss_f*Ts_avg) / (2*k + 2*(h_conv+5)*dx + 0.25*k*loss_f);
    T(N, :) = Th; 
    
    if max(abs(T - T_old), [], 'all') < 1e-6, break; end
end

%% 4. RESULTS & VISUALIZATION
sensor_cm = [4, 8, 12, 16];
s_idx = round((sensor_cm / 20) * (N - 1)) + 1;
T_final_sensors = T(s_idx, s_idx)

figure('Color', 'w');
[X, Y] = meshgrid(linspace(0, 20, N), linspace(0, 20, N));
contourf(X, Y, T, 30, 'LineColor', 'none'); 

% --- COLOR MATCHING TWEAK ---
colorbar; 
colormap jet; 
clim([40 65]); % <--- Sets the scale to match experimental range
% ----------------------------

hold on;
plot(meshgrid(sensor_cm), meshgrid(sensor_cm)', 'ko', 'MarkerSize', 10, 'LineWidth', 2);
title(['Numerical Temperature Profile ']);
xlabel('x (cm)'); ylabel('y (cm)');