
clear all, close all, clc

data = readtable('properties.xlsx','Sheet','summary');

Peval = 1; % bar
nO = 0; % moles oxygen in
nN2 = 0; % moles N2 in
ratio = 3;
egr = 100*(0:0.05:0.75);
Teval = 200:50:800; % T/ºC

[T_mesh,r_mesh] = meshgrid(Teval,egr);

for i = 1:1:length(Teval)
    for j = 1:length(egr)
    [x,gt(j,i)] = min_Gt(Peval,Teval(1,i)+273.15,ratio,nO,nN2,data,r_mesh(j,1)/100);

    X(j,i) = (1 - x(1,1)) * 100;
    H2O(j,i) = x(1,2);
    CO2(j,i) = x(1,3);
    H2(j,i) = x(1,4);
    CO(j,i) = x(1,5);
    CH4(j,i) = x(1,6);
    
    Cgr(j,i) = x(1,end-2);
    Cmw(j,i) = x(1,end-1);
    Cam(j,i) = x(1,end);
    
    C_T(j,i) = sum(x(1,end-2:end));

    end
end

set(0,'DefaultAxesFontName', 'Times New Roman')

g = figure;

tiledlayout(3,2,"TileSpacing","compact")

colormap(summer)

nexttile
surf(Teval,egr,H2)
zlabel('$\mathbf{H_{2} / \mathrm{etOH}}$', 'Interpreter', 'latex');
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex');
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.05*Teval(1), 0.95*egr(end), 0.95*max(H2(:)), '(a)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
surf(Teval, egr, CO2)
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
zlabel('$\mathbf{CO_2/\mathrm{etOH}}$', 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.05*Teval(1), 0.95*egr(end), 0.95*max(CO2(:)), '(b)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
surf(Teval, egr, CO)
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
zlabel('$\mathbf{CO/\mathrm{etOH}}$', 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.05*Teval(1), 0.95*egr(end), 0.95*max(CO(:)), '(c)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
surf(Teval, egr, CH4)
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
zlabel('$\mathbf{CH_4/\mathrm{etOH}}$', 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.05*Teval(1), 0.95*egr(end), 0.95*max(CH4(:)), '(d)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
surf(Teval, egr, H2O)
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
zlabel('$\mathbf{H_2O/\mathrm{etOH}}$', 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.05*Teval(1), 0.95*egr(end), 0.95*max(H2O(:)), '(e)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
surf(Teval, egr, C_T)
zlabel('$\mathbf{C/\mathrm{etOH}}$', 'Interpreter', 'latex')
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.05*Teval(1), 0.95*egr(end), 0.95*max(C_T(:)), '(f)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% exportgraphics(g,'EGR.png','Resolution',600)

h = figure;

tiledlayout(2,2,"TileSpacing","compact")

colormap(summer)

nexttile
pcolor(Teval, egr, C_T)
shading interp;
colorbar;
cb = colorbar;
caxis([0 0.8]);
cb.Ticks = 0:0.1:0.8;
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.1*Teval(1), 0.95*egr(end), '(a)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
pcolor(Teval, egr, Cgr)
shading interp;
colorbar;
cb = colorbar;
caxis([0 0.8]);
cb.Ticks = 0:0.1:0.8;
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.1*Teval(1), 0.95*egr(end), '(b)', 'Color', 'white', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
pcolor(Teval, egr, Cam)
shading interp;
colorbar;
cb = colorbar;
caxis([0 0.8]);
cb.Ticks = 0:0.1:0.8;
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.1*Teval(1), 0.95*egr(end), '(c)', 'Color', 'black', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

nexttile
pcolor(Teval, egr, Cmw)
shading interp;
colorbar;
cb = colorbar;
caxis([0 0.8]);
cb.Ticks = 0:0.1:0.8;
xlabel('$\mathbf{T}$ ($^\circ$C)', 'Interpreter', 'latex')
ylabel('$\mathbf{EGR\,(\%)}$', 'Interpreter', 'latex');
xlim([Teval(1) Teval(end)])
ylim([egr(1) egr(end)])
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
text(1.1*Teval(1), 0.95*egr(end), '(d)', 'Color', 'white', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% exportgraphics(h,'EGR_Cdetails.png','Resolution',600)
