function L_C = test_contour(w_so, L_B_, C_edge)
%test_contour supprime les keypoints situés sur les contours
%   ARGUMENTS:
%   	w_so: DoG digital (o=1,...,n_oct et s=0,...,n_spo+2)
%       L_B_: liste filtrée de keypoints
%       C_edge: Seuil sur le rapport entre les valeurs propres du Hessien
%   OUTPUT:
%       L_C: liste des points SIFT

%% Initialisation
L_C = [];

%% Boucle de filtrage
for i=1:size(L_B_,1)
    o = L_B_(i,1); s = L_B_(i,2); m = L_B_(i,3); n = L_B_(i,4);
    
    % Calcul du Hessien 2D
    h_11 = w_so{s+1,o+1}(m+2,n+1) + w_so{s+1,o+1}(m,n+1) - 2*w_so{s+1,o+1}(m+1,n+1);
    h_22 = w_so{s+1,o+1}(m+1,n+2) + w_so{s+1,o+1}(m+1,n) - 2*w_so{s+1,o+1}(m+1,n+1);
    h_12 = 1/4*(w_so{s+1,o+1}(m+2,n+2) - w_so{s+1,o+1}(m+2,n) - w_so{s+1,o+1}(m,n+2) + w_so{s+1,o+1}(m,n));
    h_21 = h_12;
    H_smno = [h_11,h_12;h_21,h_22];

    edgeness = trace(H_smno)^2/det(H_smno);
    if edgeness<(C_edge+1)^2/C_edge;
        L_C = [L_C; L_B_(i,:)];
    end
end

end

