
% Objective function

function G = phase_equilibrium(nj,Peval,T,u_0,data)

P0 = 1; % bar

R = 1.987; % ca/mol/K

Enj = sum(nj(1:end-3));

y = nj/Enj; y = y';

phi = PR_EOS_vdW(T,Peval,data,y(1:length(nj)-3));

P = phi.*(Peval/P0); % bar/bar

G = 0;
for i = 1:length(nj)-3
    G = G + nj(i) * ( u_0(i) + R * T * ( log(P(i)) + log(y(i)) ) );
end

for j = length(nj)-2:length(nj)
    G = G + nj(j) * u_0(j);
end

end