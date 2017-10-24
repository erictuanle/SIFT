function [partial_m, partial_n] = gradient_2D(v_so, M, N, n_oct, n_spo)
%gradient_2D calcule le gradient 2D de chaque image du scale space
%   ARGUMENTS:
%   	v_so: scale space digital (o=1,...,n_oct et s=0,...,n_spo+2)
%   OUTPUT:
%       patial_m: gradient du scale-space selon x (o=1,...,n_oct et s=1,...,n_spo)
%       patial_n: gradient du scale-space selon y (o=1,...,n_oct et s=1,...,n_spo)
for o=1:n_oct
    Mo = floor(M/2^(o-1));
    No = floor(N/2^(o-1));
    for s=1:n_spo
        temp_m = zeros(Mo-2,No-2);
        temp_n = zeros(Mo-2,No-2);
        for m=1:Mo-2
            for n=1:No-2
                temp_m(m,n) = 1/2*(v_so{s+1,o+1}(m+2,n+1)-v_so{s+1,o+1}(m,n+1));
                temp_n(m,n) = 1/2*(v_so{s+1,o+1}(m+1,n+2)-v_so{s+1,o+1}(m+1,n));
            end
        end
        partial_m{s+1,o+1} = temp_m;
        partial_n{s+1,o+1} = temp_n;
    end
end

end

