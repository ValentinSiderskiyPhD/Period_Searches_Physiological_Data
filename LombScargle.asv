clc
clear
close all
%% Problem 1b

x = linspace(0,10,123);  %array-like input values to evaluate the array
A= 2;  %Amplitude
p = pi;  %period
phi = 0; %phase
noise = 0; %noise amplitude
y = gen_periodic_data(x,p,A,phi,noise);

figure('Name','Problem 1b: Test gen_periodic_data function')
plot(x,y,'*')


%% Problem 1d:Test phase_plot
pp_fig = phase_plot(x,y, pi);
set(gcf,'Name','Problem 1d: Test phase_plot')

%% Problem 1d+:Test phase_plot with error bars
[m,n] = size(x);
noise = 3;
ppp_fig = phase_plot(x,y, pi,wgn(m,n,noise,'linear'));
set(gcf,'Name','Problem 1d+: Test phase_plot with error bars')


%% Problem 2a
x = 0:10;
P = 2;
A = 2;
noise = 0;
y = gen_periodic_data(x,P,A,phi,noise);
x_signal = linspace(0,10,1234);
y_signal = gen_periodic_data(x_signal,P,A,phi,noise);
figure('Name','Problem 2a: samplint at the Nyquist frequency')
plot(x,y,'*')
hold on
plot(x_signal,y_signal)
hold off

%% Problem 2b
x = 0:10;
f = 0.7;
P = 1/(f);
y = gen_periodic_data(x,P,A,phi,noise);
x_signal = linspace(0,10,1234);
y_signal = gen_periodic_data(x_signal,P,A,phi,noise);
figure('Name','Problem 2b: sampling at blow the Nyquist frequency')
plot(x,y,'*')
hold on
plot(x_signal,y_signal)


%% Problem 2c
x = 0:10;
f = 2.7;
P = 1/(f);

y_high = gen_periodic_data(x,P,A,phi,noise);
y_high_signal = gen_periodic_data(x_signal,P,A,phi,noise);
figure('Name','Problem 2c: sampling much blow the Nyquist frequency')

plot(x,y,'*')
hold on
plot(x_signal,y_signal)
plot(x,y_high,'*')
plot(x_signal,y_high_signal)
hold off


%% Problem 3c+ test the min_chi2 function
A = 7.4;
P = 5.25;
phi = pi/4;
noise = 0.8;
x = 10*rand(100,1);
y = gen_periodic_data(x,P,A,phi,noise);
y_unc = ones(size(x))*sqrt(0.8);

% Note in the lecture P (period) and f (frequency) get swapped around in the code. Confusing.
P_grid = linspace(1/peak2peak(x),10,50);
theta = [0,0];

theta_values = zeros(2,1,length(P_grid));
chi2_values = zeros(size(P_grid));

figure('Name', 'Check that min chi2 is actually a minimum, and A and phi estimates are good')
hold on 
for i = 1:length(P_grid)
    [theta_values(:,:,i),chi2_values(i)] =  min_chi2(theta, y, y_unc, x, 1./P_grid(i));
end

[M,I] = min(chi2_values);
A_opt = theta_values(1,1,I);
phi_opt= theta_values(2,1,I);

plot(P_grid,chi2_values,'*');

title( ['A = ', num2str(A),', A_{opt} = ', num2str(A_opt), ' phi = ', num2str(phi), ' phi_{opt} = ', num2str(phi_opt)] )
theta_values(2,1,I);

%% Problem 3d

A = 7.4;
P = 5.25;
phi = 0;
noise = 0.8;

x_signal = linspace(0,10,1234);
y_signal = gen_periodic_data(x_signal,P,A,phi,noise);

x = 10*rand(100,1);
y = gen_periodic_data(x,P,A,phi,noise);
y_unc = ones(size(x))*sqrt(0.8);

f_grid = linspace(1/peak2peak(x),10,50);

psd_ls = ls_periodogram(y, y_unc, x, f_grid);
figure('Name','ls_perodogram test')
plot(1./f_grid, psd_ls)
xlabel('period')
ylabel('P')

hold on 

f_grid_fine = linspace(1/peak2peak(x),10,1000);
psd_ls_fine = ls_periodogram(y, y_unc, x, f_grid_fine);
plot(1./f_grid_fine, psd_ls_fine)
xlabel('period')
ylabel('P')

%% Problem 3f

[~,I] = max(psd_ls_fine);
f_opt = f_grid_fine(I);
phase_plot(x,y, 1/f_opt, y_unc);

%% Proglem 3g

f_min = 1/(10*365);
f_max = 1/(1/24);
delta_f = f_min/5;

f_grid = f_min:delta_f:f_max;
size(f_grid)



%% Problem 1a
function periodic_data = gen_periodic_data(x,period,amplitude,phase,noise)

    [m,n] = size(x);
    dy = wgn(m,n,noise,'linear');
    y = amplitude *sin(2*pi*x/period - phase);
    periodic_data = y + dy;

end

%% Problem 1c

% Takes x,y, P as input to create a phase-folded light curve (i.e. plot the
% data at their respective phase values given the period P).

% Include an optonal arument y_unc.
function fig_h = phase_plot(varargin)

    if( nargin ==3 || nargin == 4)
        x = varargin{1};
        y = varargin{2};
        period = varargin{3};
    else
        error('Unexpected inputs')    
    end


    switch nargin
        case 4  % x,y,period,y_unc
            y_unc = varargin{4};
        case 3  % x,y,period
            y_unc = zeros(1,length(x));
        otherwise
            error('Unexpected inputs')    
    end

    phases = mod((x/period),1);
    
    [~,plot_order] = sort(phases);
    figure
    fig_h = errorbar(phases(plot_order),y(plot_order),y_unc,'*');
    set(gca, 'XLimSpec', 'Tight');

end

%% Proglem 3a  

function  chi2 = getChi2(theta, y, y_unc, x, f)

    a= theta(1);
    phi = theta(2);
    chi2 = sum(((y- a*sin(2*pi*x*f - phi)).^2)./(y_unc.^2));

end 

%% Problem 3b  wite a function to minimize chi2
function [theta_opt,chi2_opt] =  min_chi2(theta, y, y_unc, x, f)

    fun = @(theta)getChi2(theta,y, y_unc, x, f);
    [theta_opt,chi2_opt] = fminsearch(fun,theta); 

end

%% Problem 3c

function psd = ls_periodogram(y, y_unc, x, f_grid)

    psd = zeros(size(f_grid));
    chi2_0 = sum(((y-mean(y))./y_unc).^2);

    for i = 1:length(f_grid)
        [~,chi2_opt] = min_chi2( [0,0],y,y_unc, x, f_grid(i));
        psd(i) = 0.5*(chi2_0 - chi2_opt);
    end

end