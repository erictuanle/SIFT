function L_A_ = test_conservatif(w_so, C_DoG, L_A)
%test_conservatif supprime les keypoints non contrastés
%   ARGUMENTS:
%   	w_so: DoG digital (o=1,...,n_oct et s=0,...,n_spo+2)
%       C_DoG: Seuil de filtrage
%       L_A: liste non filtrée des extremas du DoG
%   OUTPUT:
%       L_A_: liste filtrée des extremas du DoG

%% Initialisation
L_A_ = [];

%% Boucle de filtrage
for i=1:size(L_A,1)
    o = L_A(i,1); s = L_A(i,2);
    m = L_A(i,3); n = L_A(i,4);
    if w_so{s+1,o+1}(m+1,n+1)>=0.8*C_DoG
        L_A_ = [L_A_; L_A(i,:)];
    end
end

end

