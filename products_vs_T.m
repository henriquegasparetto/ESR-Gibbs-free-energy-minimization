
clear all, close all, clc

data = readtable('properties.xlsx','Sheet','summary');

Peval = 1; % bar
nO = 0; % moles oxygen in
nN2 = 0; % moles N2 in
ratio = 3; % H2O/ethanol/mol/mol
egr = 0;
Teval = 200:25:800; % T/ºC

[T_mesh,r_mesh] = meshgrid(Teval,ratio);

for i = 1:1:length(Teval)
    [x,gt(1,i)] = min_Gt(Peval,Teval(1,i)+273.15,ratio,nO,nN2,data,egr);
    
    X(1,i) = (1 - x(1,1)) * 100;
    H2O(1,i) = x(1,2);
    CO2(1,i) = x(1,3);
    H2(1,i) = x(1,4);
    CO(1,i) = x(1,5);
    CH4(1,i) = x(1,6);
    others(1,i) = sum(x(1,7:end-3));
    C_T(1,i) = sum(x(1,end-2:end));
    
end

resultsTable = table(Teval', X', H2O', CO2', H2', CO', CH4', others', C_T', gt', ...
    'VariableNames', {'T', 'X', 'H2O', 'CO2', 'H2', 'CO', 'CH4', 'Others', 'C_T', 'Gt'});

% writetable(resultsTable, 'ESR_T_ratio5.6_1bar.csv');

set(0,'DefaultAxesFontName', 'Times New Roman')

colororder(linspecer(6))

plot(Teval, H2, '-', 'LineWidth', 2);
hold on;
plot(Teval, CO, '-', 'LineWidth', 2);
plot(Teval, CO2, '-', 'LineWidth', 2);
plot(Teval, CH4, '-', 'LineWidth', 2);
plot(Teval, C_T, '-', 'LineWidth', 2);
plot(Teval, H2O, '-', 'LineWidth', 2);
hold off;

legend({'H$_{2}$', 'CO', 'CO$_{2}$', 'CH$_{4}$', 'C', 'H$_{2}$O'},'Location','best','Interpreter','latex')
legend('box','off')
ylabel('$\mathrm{Yield} \, (\mathrm{mol/mol \, etOH})$','Interpreter','latex');
xlabel('$\mathbf{T}$ ($^\circ$C)','Interpreter','latex')
xlim([Teval(1) Teval(end)])
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'LineWidth', 1.5, ...
    'TickLength', [0.02 0.02], 'Box', 'on');

% exportgraphics(gca,'ESR_T_ratio5.6_1bar.png','Resolution',600)
