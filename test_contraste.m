function L_B_ = test_contraste(L_B, C_DoG)
%test_contraste permet de retirer les keypoints peu contrastés
%   ARGUMENTS:
%   	L_B: liste non filtrée de keypoints
%       C_DoG: seuil de filtrage
%   OUTPUT:
%       L_B_: liste filtrée de keypoints

%% Initialisation
L_B_ = [];

%% Boucle de filtrage
for i=1:size(L_B,1)
    omega = L_B(i,8);
    if abs(omega)>=C_DoG
        L_B_ = [L_B_; L_B(i,:)];
    end
end

end