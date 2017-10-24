function [A_dict, b_dict, copie_dict] = SATS(matching, t1, t2)
% test_RANSAC permet de vérifier les matchs par invariance par rotation
%   et changement d'échelle
% ARGUMENTS:
%   matching: liste des matchs
%   t1: seuil 1
%   t2: seuil 2
% OUTPUT:
%   A: dictionnaire des matrices des transformations affines inter-patchs
%   b: dictionnaire des translations inter-patchs
%   copie: dictionnaire des zones copiées

cpt = 1;
for i=1:size(matching,1)
    if mod(i,100)==0
        disp([num2str(floor(i/size(matching,1)*100)),'%'])
    end
    
    H = matching(i,:)';
    for j=1:size(matching,1)
        if i~=j && norm(matching(i,1:2)-matching(j,1:2))<t1 && norm(matching(i,3:4)-matching(j,3:4))<t1
            H = [H,matching(j,:)'];
        end
    end
    
    % Pour déterminer A, il faut s'assurer que l'on ait bien au moins trois
    %   matchs
    if size(H,2)>=3
        % Calcul de A et b
        try
            aff_2D = estimateGeometricTransform(H(1:2,:)',H(3:4,:)','affine');
            mat = aff_2D.T;
            A = mat(1:2,1:2);
            b = mat(3,1:2);

            for j=1:size(matching,1)
                distances = [];
                for k=1:size(H,2)
                    distances = [distances, norm(matching(j,1:2)-H(1:2,k)')];
                end
                distances = (distances<t1) .* (distances>0);
                
                if i~=j && sum(distances)>=1
                    p_i = matching(j,1:2);
                    q_i = p_i*A+b;
                    if norm(matching(j,3:4)-q_i)<t1
                        H = [H,[matching(j,1:2)';matching(j,3:4)']];
                        H = unique(H.','rows').';
                        if mod(size(H,2),10)==0
                            aff_2D = estimateGeometricTransform(H(1:2,:)',H(3:4,:)','affine');
                            mat = aff_2D.T;
                            A = mat(1:2,1:2);
                            b = mat(3,1:2);
                        end
                    end
                end
            end

            % Sauvegarde de A, b et des blocs copiés de H
            if size(H,2)>=t2
                A_dict{cpt} = A;
                b_dict{cpt} = b;
                copie_dict{cpt} = H;         
                cpt = cpt+1;
            end
        catch
            cpt = cpt;
        end
    end
end

if  ~exist('A_dict')
    A_dict = [];
    b_dict = [];
    copie_dict = [];
end

end