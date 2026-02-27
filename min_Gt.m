
function [x,Gt] = min_Gt(Peval,T,ratio,nO,nN2,data,egr)

T0 = 298; % K

in = readmatrix('properties.xlsx','Sheet','in'); 

out = readmatrix('properties.xlsx','Sheet','out'); out = out(:,2:end)';

ratio = ratio * (1 - egr);

aux = (1 + ratio) * egr; % total of EGR moles

in(1,:) = (1 - egr) * in(1,:); % ethanol
in(2,:) = (aux * 0.1921) + (ratio * in(2,:)); % H2O
in(3,:) = (aux * 0.0650) + in(3,:); % CO2
in(5,:) = (aux * 0.0285) + in(5,:); % CO
in(12,:) = (aux * 0.0922) + nO * in(12,:); % O2
in(15,:) = (aux * 0.6222) + nN2 * in(15,:); % N2

in = sum(in); in = in(:,2:end)';

u_0 = (T / T0) * data.Gf ...
     + (1 - (T / T0)) * data.Hf ...
     - data.CpA * (T * log(T / T0) - T + T0) ...
     - (data.CpB/2) * (T^2 - 2 * T * (T0^1) + 1 * (T0^1)) ...
     - (data.CpC/6) * ( T^3 - 3 * T * (T0^2) + 2 * (T0^3)) ...
     - (data.CpD/12) * ( T^4 - 4 * T * (T0^3) + 3 * (T0^4)) ...
     - (data.CpE/20) * ( T^5 - 5 * T * (T0^4) + 4 * (T0^5));
mwcnt = u_0(end) + (8250 - 11.72 * T)/4.184;
amorphous = (- (5.8238e-12) * (T^4) + (1.9769e-8) * (T^3) ...
     - (2.7622e-5) * (T^2) + (9.8415e-3) * T + 14.895)/4.184;

NN = length(u_0);
 
u_0(NN+1,1) = mwcnt;
u_0(NN+2,1) = amorphous;

out(:,NN+1) = out(:,NN);
out(:,NN+2) = out(:,NN);

u_0 = u_0'; % cal/mol

LB = zeros(1,length(u_0)); % no mole numbers less than zero

options = optimset('Algorithm','sqp');
[x,Gt] = fmincon(@phase_equilibrium,rand(1,length(u_0)),[],[],out,in,LB,[],[],options,Peval,T,u_0,data);

end