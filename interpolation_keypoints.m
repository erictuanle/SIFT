function L_B = interpolation_keypoints(L_A_, w_so, delta_min, sigma_min, n_spo, M, N)
%interpolation_keypoints interpole les keypoints dans l'esapce 3D
%   ARGUMENTS:
%       L_A_: : liste filtrée des extremas du DoG (discret)
%   	w_so: DoG digital (o=1,...,n_oct et s=0,...,n_spo+2)
%   OUTPUT:
%       L_B: liste des extremas du DoG (interpolé)

%% Initialisation
L_B = [];

%% Boucle d'interpolation
for i=1:size(L_A_,1)
    o_e = L_A_(i,1); s_e = L_A_(i,2); m_e = L_A_(i,3); n_e = L_A_(i,4);
    Mo_e = floor(M/2^(o_e-1));
    No_e = floor(N/2^(o_e-1));
    s = s_e; m = m_e; n = n_e;
    cpt = 0;
    alpha_star = inf(3,1);
    while max(abs(alpha_star(1)),max(abs(alpha_star(2)),abs(alpha_star(3))))>=0.6 && cpt<5
        [alpha_star, omega] = interpolation_quadratique(w_so,o_e,s,m,n);
        
        delta_oe = delta_min*2^(o_e-1);
        sigma = delta_oe/delta_min*sigma_min*2^((alpha_star(1)+s)/n_spo);
        x = delta_oe*(alpha_star(2)+m); y = delta_oe*(alpha_star(3)+n);
        
        s = round(s+alpha_star(1)); m = round(m+alpha_star(2)); n = round(n+alpha_star(3));
        cpt = cpt+1;
        if ~(s>=1 && s<=n_spo && m>=1 && m<=Mo_e-2 && n>=1 && n<=No_e-2)
            cpt = 5;
        end
    end
    
    if max(abs(alpha_star(1)),max(abs(alpha_star(2)),abs(alpha_star(3))))<0.6
        if s>=1 && s<=n_spo && m>=1 && m<=Mo_e-2 && n>=1 && n<=No_e-2
            L_B = [L_B;o_e,s,m,n,sigma,x,y,omega];
        end
    end
end

end