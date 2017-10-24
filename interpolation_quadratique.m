function [alpha_star, omega] = interpolation_quadratique(w_so, o, s, m, n)
%interpolation_quadratique effectue l'interpolation quadratique des
%       keypoints
%   ARGUMENTS:
%   	w_so: DoG digital (o=1,...,n_oct et s=0,...,n_spo+2)
%       (o,s,m,n): coordonnées de l'extremum 3D du DoG
%   OUTPUT:
%       alpha_star: offset du centre de l'extremum 3D interpolé
%       omega: valeur de l'extremum interpolé

%% Calcul du gradient et du Hessien
g_smno = [1/2*(w_so{s+2,o+1}(m+1,n+1)-w_so{s,o+1}(m+1,n+1));...
    1/2*(w_so{s+1,o+1}(m+2,n+1)-w_so{s+1,o+1}(m,n+1));...
    1/2*(w_so{s+1,o+1}(m+1,n+2)-w_so{s+1,o+1}(m+1,n))];

h_11 = w_so{s+2,o+1}(m+1,n+1) + w_so{s,o+1}(m+1,n+1) - 2*w_so{s+1,o+1}(m+1,n+1);
h_12 = 1/4*(w_so{s+2,o+1}(m+2,n+1) - w_so{s+2,o+1}(m,n+1) - w_so{s,o+1}(m+2,n+1) + w_so{s,o+1}(m,n+1));
h_22 = w_so{s+1,o+1}(m+2,n+1) + w_so{s+1,o+1}(m,n+1) - 2*w_so{s+1,o+1}(m+1,n+1);
h_13 = 1/4*(w_so{s+2,o+1}(m+1,n+2) - w_so{s+2,o+1}(m+1,n) - w_so{s,o+1}(m+1,n+2) + w_so{s,o+1}(m+1,n));
h_33 = w_so{s+1,o+1}(m+1,n+2) + w_so{s+1,o+1}(m+1,n) - 2*w_so{s+1,o+1}(m+1,n+1);
h_23 = 1/4*(w_so{s+1,o+1}(m+2,n+2) - w_so{s+1,o+1}(m+2,n) - w_so{s+1,o+1}(m,n+2) + w_so{s+1,o+1}(m,n));
H_smno = [h_11,h_12,h_13;h_12,h_22,h_23;h_13,h_23,h_33];

%% Calcul de alpha_star
alpha_star = -inv(H_smno)*g_smno;
omega = w_so{s+1,o+1}(m+1,n+1) - 1/2*g_smno'*inv(H_smno)*g_smno;
end