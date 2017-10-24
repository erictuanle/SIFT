function u_ = interpolation_bilineaire(u, delta_)
%interpolation_bilineaire effectue l'interpolation bilinéaire d'une image
%   ARGUMENTS:
%   	u: image digitale de taille MxN
%       delta_: distance inter-pixel de l'image de sortie
%   OUTPUT:
%       u_: image digitale de taille M'xN' avec M'=floor(M/delta_) et N'=floor(N/delta_)

%% Paramètres
[M,N] = size(u);
M_ = floor(M/delta_);
N_ = floor(N/delta_);
u_ = zeros(M_,N_);

%% Boucle de calcul
for m_=0:M_-1
    for n_=0:N_-1
        x = delta_*m_;
        y = delta_*n_;
        u_(m_+1,n_+1) = (x-floor(x))*(y-floor(y))*...
            u(min(mod(1+floor(x),2*M),mod(2*M-2-floor(x),2*M))+1,min(mod(1+floor(y),2*N),mod(2*N-2-floor(y),2*N))+1)...
            +(1+floor(x)-x)*(y-floor(y))*...
            u(min(mod(floor(x),2*M),mod(2*M-1-floor(x),2*M))+1,min(mod(1+floor(y),2*N),mod(2*N-2-floor(y),2*N))+1)...
            +(x-floor(x))*(1+floor(y)-y)*...
            u(min(mod(floor(x)+1,2*M),mod(2*M-2-floor(x),2*M))+1,min(mod(floor(y),2*N),mod(2*N-1-floor(y),2*N))+1)...
            +(1+floor(x)-x)*(1+floor(y)-y)*...
            u(min(mod(floor(x),2*M),mod(2*M-1-floor(x),2*M))+1,min(mod(floor(y),2*N),mod(2*N-1-floor(y),2*N))+1);
    end
end

end

