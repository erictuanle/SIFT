function v_so = scale_space_gaussien(u_in, n_oct, n_spo, sigma_min, delta_min, sigma_in)
%gaussian_scale_space Calcul du scale space digital Gaussien
%   ARGUMENTS:
%       u_in: image digitale en entrée de taille MxN
%   OUTPUT:
%       v_so: scale space digital (o=1,...,n_oct et s=0,...,n_spo+2)

%% Paramètres
[M,N] = size(u_in);

%% Calcul de la première octave
u_ = interpolation_bilineaire(u_in, delta_min);
sigma = 1/delta_min*sqrt(sigma_min^2-sigma_in^2);
v_so{1,2} = imgaussfilt(u_,sigma);
for s=1:n_spo+2
    sigma = sigma_min/delta_min*sqrt(2^(2*s/n_spo)-2^(2*(s-1)/n_spo));
    v_so{s+1,2} = imgaussfilt(v_so{s,2},sigma);
end

%% Calcul des octaves suivantes
for o=2:n_oct
    Mo = floor(M/2^(o-1));
    No = floor(N/2^(o-1));
    
    v_temp = zeros(Mo,No);
    for m=0:Mo-1
       for n=0:No-1
           v_temp(m+1,n+1) = v_so{n_spo+1,o}(2*m+1,2*n+1);
       end
    end
    v_so{1,o+1} = v_temp;
    
    for s=1:n_spo+2
       sigma = sigma_min/delta_min*sqrt(2^(2*s/n_spo)-2^(2*(s-1)/n_spo));
       v_so{s+1,o+1} = imgaussfilt(v_so{s,o+1},sigma);
    end
end

end