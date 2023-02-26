%Plotting the Test Data
clc
clear
close all
load('VVA012_TILT0_03_STACK.mat')

start_stack = 2; %ignore the first stack

ind0 = length(stack{1});

%MCA cm
for i = 1:length(ind0)
    
    figure('Name','MCA')
    
    subplot(2,1,1)
    ylabel('tilt (deg)')
    hold on
    for j = start_stack:length(stack{i})
        plot(stack{i}{j}.time-stack{i}{j}.time(1),stack{i}{j}.Tilt,'DisplayName',['Osci#', num2str(j)]);
    end    
    legend
    set(gca, 'XLimSpec', 'Tight');
    
    subplot(2,1,2)
    ylabel('MCA velocity (cm/s)')
    hold on
    
    for j = start_stack:length(stack{i})
        
        a = plot(stack{i}{j}.time-stack{i}{j}.time(1),stack{i}{j}.MBVcm,'DisplayName',['Oscillation #', num2str(j)]);
        legend
        plot(tRRi(stack{i}{j}.tRRi_idx)-stack{i}{j}.time(1),MBVcm( stack{i}{j}.tRRi_idx),'x','Color', a.Color, 'HandleVisibility','off')
    end
       
    set(gca, 'XLimSpec', 'Tight');
    
end
