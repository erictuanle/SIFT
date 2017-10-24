function L_E = descripteurs_keypoints(L_D, partial_m, partial_n, n_hist, n_ori, lambda_descr, delta_min)
%orientation_keypoints calcule les descripteurs associés à chaque keypoint
%   ARGUMENTS:
%   	L_D: liste des points SIFT orientés
%       patial_m: gradient du scale-space selon x (o=1,...,n_oct et s=1,...,n_spo)
%       patial_n: gradient du scale-space selon y (o=1,...,n_oct et s=1,...,n_spo)
%       n_hist: nombre d'histogrammes
%       n_ori: nombre de bins dans les histogrammes
%       lambda_descr: paramètres de la fenêtre gaussienne
%   OUTPUT:
%       L_E: liste des keypoints avec descripteurs
L_E = [];
for elt=1:size(L_D,1)
    o_key = L_D(elt,1); s_key = L_D(elt,2); x_key = L_D(elt,3); y_key = L_D(elt,4);
    sigma_key = L_D(elt,5); theta_key = L_D(elt,7);

    % Vérification si le point est bien assez loi du bord
    if round((x_key-sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1)))>=1 &&...
            round((x_key+sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1)))<=size(partial_m{s_key+1,o_key+1},1) &&...
            round((y_key-sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1)))>=1&&...
            round((y_key+sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1)))<=size(partial_m{s_key+1,o_key+1},2)
        
        % Initialisation de l'histogramme
        h = zeros(n_hist,n_hist,n_ori);
        
        % Echantillons du patch normalisaé P^descr
        m_list = round((x_key-sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1))):round((x_key+sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1)));
        n_list = round((y_key-sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1))):round((y_key+sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/(delta_min*2^(o_key-1)));
        for m=m_list
            for n=n_list
                x_mn = 1/sigma_key*((m*delta_min*2^(o_key-1)-x_key)*cos(theta_key)+(n*delta_min*2^(o_key-1)-y_key)*sin(theta_key));
                y_mn = 1/sigma_key*(-(m*delta_min*2^(o_key-1)-x_key)*sin(theta_key)+(n*delta_min*2^(o_key-1)-y_key)*cos(theta_key));
                
                if max(abs(x_mn),abs(y_mn))<lambda_descr*(n_hist+1)/n_hist
                    theta_mn = mod(atan2(partial_n{s_key+1,o_key+1}(m,n),partial_m{s_key+1,o_key+1}(m,n))-theta_key,2*pi);
                    c_mn = exp(-sqrt((m*delta_min*2^(o_key-1)-x_key)^2+(n*delta_min*2^(o_key-1)-y_key)^2)^2/(2*(lambda_descr*sigma_key)^2))*sqrt((partial_m{s_key+1,o_key+1}(m,n))^2+(partial_n{s_key+1,o_key+1}(m,n))^2);
                    for i=1:n_hist
                        for j=1:n_hist
                            x_i = (i-(1+n_hist)/2)*2*lambda_descr/n_hist;
                            y_i = (j-(1+n_hist)/2)*2*lambda_descr/n_hist;
                            if abs(x_i-x_mn)<=2*lambda_descr/n_hist && abs(y_i-y_mn)<=2*lambda_descr/n_hist
                                for k=1:n_ori
                                   theta_k = 2*pi*(k-1)/n_ori;
                                   if abs(mod(theta_k-theta_mn,2*pi))<2*pi/n_ori
                                       h(i,j,k) = h(i,j,k) + (1-n_hist/(2*lambda_descr)*abs(x_mn-x_i))*(1-n_hist/(2*lambda_descr)*abs(y_mn-y_i))*...
                                           (1-n_ori/(2*pi)*abs(mod(theta_mn-theta_k,2*pi)))*c_mn;
                                   end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % Construction du vecteur de features
        f = zeros(1,n_hist*n_hist*n_ori);
        for i=1:n_hist
            for j=1:n_hist
                for k=1:n_ori
                    f((i-1)*n_hist*n_ori+(j-1)*n_ori+k)=h(i,j,k);
                end
            end
        end
        
        % Quantification du vecteur de features
        for l=1:n_hist*n_hist*n_ori
            f(l) = min(f(l),0.2*norm(f));
            f(l) = min(floor(512*f(l)/norm(f)),255);
        end
        
        L_E = [L_E;x_key,y_key,sigma_key,theta_key,f];
    end
end

end