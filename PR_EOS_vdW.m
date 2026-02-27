
% Peng-Robinson equation of state with van der Walls mixing rule

function phi = PR_EOS_vdW(T,P,data,y)

data = data(1:end-1,:); % remove C

R_PR = 8.314472; % J/mol/K

m = zeros(length(data.Tc_K),1);
for i = 1:length(data.Tc_K)
    if data.w(i) <= 0.491
        m(i,1) = 0.37464 + 1.54226 * data.w(i) - 0.26992 * (data.w(i)^2);
    else
        m(i,1) = 0.379642 + 1.48503 * data.w(i) - 0.164423 * (data.w(i)^2) + 0.016666 * (data.w(i)^3);
    end
end

a = (0.457235529 * (R_PR^2) .* (data.Tc_K.^2) ./ (data.Pc_bar * 100000)) .* (1 + m .* (1 - sqrt(T ./ data.Tc_K))).^2;

b = 0.0777960739 * R_PR .* data.Tc_K ./ (data.Pc_bar * 100000);

bm = sum(y .* b);

% rough approximation for kij using Gao et al. correlation for BIP (predicted for light HC)
% necessary for mixing rules
for i = 1:length(data.Tc_K)
    for j = 1:length(data.Tc_K)
        Zc(j,i) = (data.Zc(i) + data.Zc(j))/2;
        
        k(j,i) = 1 - ((2 * sqrt(data.Tc_K(i) * data.Tc_K(j)) / (data.Tc_K(i) + data.Tc_K(j)))^Zc(j,i));
    end
end

am = 0;
for i = 1:length(a)
    for j = 1:length(a)
        am = am + y(i) * y(j) * sqrt(a(i) * a(j)) * (1 - k(i,j));
    end
end

% for a mixture there is only one Z so am and bm must be used

A = am .* (P * 100000) ./ (R_PR^2 * T^2);

B = bm .* (P * 100000) ./ (R_PR * T);

% Z^3 - (1 - B) * (Z^2) + (A -3*(B^2) -2*B)*Z - (A*B -B^2 -B^3) = 0

c2 = -(1 - B);
c1 = A - 3 * B.^2 - 2*B;
c0 = -(A .* B - B.^2 - B.^3);

Z = zeros(length(A),1);
for i = 1:length(A)
    coeffs = [1, c2(i), c1(i), c0(i)];
    Z_roots = roots(coeffs);
    Z_roots = Z_roots(imag(Z_roots) == 0);
    Z(i,1) = max(Z_roots);
end

% aux variable to compute fugacity coefficient
aux = zeros(length(a),1);
for i = 1:length(a)
    for j = 1:length(a)
        aux(i) = aux(i) + y(j) * sqrt(a(i) * a(j)) * (1 - k(i,j));
    end
end
aux = 2 * aux;

phi = exp(((b/bm) .* (Z - 1)) - log(Z - B) + (A ./ (2 * sqrt(2) .* B)) .* ...
    (aux/am - b/bm) .* ...
    log((Z + (1 + sqrt(2)) .* B) ./ (Z + (1 - sqrt(2)) .* B)));

end
