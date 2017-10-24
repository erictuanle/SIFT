function w_so = DoG(v_so, n_oct, n_spo)
%DoG calcule la Difference of Gaussians scale-space
%   ARGUMENTS:
%   	v_so: scale space digital (o=1,...,n_oct et s=0,...,n_spo+2)
%   OUTPUT:
%       w_so: DoG digital (o=1,...,n_oct et s=0,...,n_spo+2)

%% Boucle de calcul
for o=1:n_oct
    for s=0:n_spo+1
        w_so{s+1,o+1} = v_so{s+2,o+1} - v_so{s+1,o+1};
    end
end

end

