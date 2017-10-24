function L_D = orientation_keypoints(L_C, partial_m, partial_n, lambda_ori, n_bins, t, delta_min)
%orientation_keypoints détermine l'orientation du patch centré autour du
%       keypoint
%   ARGUMENTS:
%   	L_C: liste des points SIFT
%       patial_m: gradient du scale-space selon x (o=1,...,n_oct et s=1,...,n_spo)
%       patial_n: gradient du scale-space selon y (o=1,...,n_oct et s=1,...,n_spo)
%       lambda_ori: paramètre pour le patch
%       n_bins: nombre de bins pour l'histogramme des orientations h
%       t: seuil pour l'orientation secondaire
%   OUTPUT:
%       L_D: liste des keypoints orientés

L_D = [];
for i=1:size(L_C,1)
    o_key = L_C(i,1); s_key = L_C(i,2); x_key = L_C(i,6); y_key = L_C(i,7);
    sigma_key = L_C(i,5); omega = L_C(i,8);

    % Vérification si le point est bien assez loi du bord
    if round((x_key-3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1)))>=1 &&...
            round((x_key+3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1)))<=size(partial_m{s_key+1,o_key+1},1) &&...
            round((y_key-3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1)))>=1&&...
            round((y_key+3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1)))<=size(partial_m{s_key+1,o_key+1},2)
    
        % Initialisation de l'histogramme
        h = zeros(n_bins,1);
        
        % Echantillons du patch normalisaé P^ori
        m_list = round((x_key-3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1))):round((x_key+3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1)));
        n_list = round((y_key-3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1))):round((y_key+3*lambda_ori*sigma_key)/(delta_min*2^(o_key-1)));
        for m=m_list
            for n=n_list
                c = exp(-sqrt((m*delta_min*2^(o_key-1)-x_key)^2+(n*delta_min*2^(o_key-1)-y_key)^2)^2/(2*(lambda_ori*sigma_key)^2))*sqrt((partial_m{s_key+1,o_key+1}(m,n))^2+(partial_n{s_key+1,o_key+1}(m,n))^2);
                bin = round(n_bins/(2*pi)*mod(atan2(partial_n{s_key+1,o_key+1}(m,n),partial_m{s_key+1,o_key+1}(m,n)),2*pi));
                if bin == 0
                    bin = n_bins;
                end
                h(bin) = h(bin) + c;
            end
        end
        
        % Smoothing de h
        for cpt=1:6
            h = cconv(h,1/3*[1;1;1],n_bins);
        end
        
        % Extraction des orientations de références
        for k=1:n_bins
            theta_k = 2*pi*(k-1)/n_bins;
            if k==1
                test = h(k)>h(n_bins) && h(k)>h(k+1) && h(k)>t*max(h);
            elseif k==n_bins
                test = h(k)>h(k-1) && h(k)>h(1) && h(k)>t*max(h);
            else
                test = h(k)>h(k-1) && h(k)>h(k+1) && h(k)>t*max(h);
            end
            if test==1 && k==1
                theta_key = theta_k + pi/n_bins*(h(n_bins)-h(k+1))/(h(n_bins)-2*h(k)+h(k+1));
                L_D = [L_D; o_key, s_key, x_key, y_key, sigma_key, omega, theta_key] ;
            elseif test==1 && k==n_bins
                theta_key = theta_k + pi/n_bins*(h(k-1)-h(1))/(h(k-1)-2*h(k)+h(1));
                L_D =  [L_D; o_key, s_key, x_key, y_key, sigma_key, omega, theta_key] ;
            elseif test==1 && k~=1 && k~=n_bins
                theta_key = theta_k + pi/n_bins*(h(k-1)-h(k+1))/(h(k-1)-2*h(k)+h(k+1));
                L_D =  [L_D; o_key, s_key, x_key, y_key, sigma_key, omega, theta_key] ;
            end
        end
    end
end

end