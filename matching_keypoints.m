function M = matching_keypoints_v2(L_E, C_matchr, delta_min)
%matching_keypoints_v2 match les keypoints d'une même image
%   ARGUMENTS:
%   	L_E: liste des descripteurs SIFT
%		C_matchr: seuil relatif de match
%   OUTPUT:
%       M: correspondance
M = [];
for elt=1:size(L_E,1)
    distances = [];
    
    cpt = 1;
    for elt_comp=1:size(L_E,1)
        if abs(L_E(elt,1)/delta_min-L_E(elt_comp,1)/delta_min)>3 && abs(L_E(elt,2)/delta_min-L_E(elt_comp,2)/delta_min)>3
            vec = [cpt ; norm(L_E(elt,5:end)-L_E(elt_comp,5:end))];
            distances = [distances, vec];
        end
        cpt = cpt+1;
    end
    
    [val,idx] = sort(distances(2,:),'ascend');
    if val(1)<C_matchr*val(2)
        idx = idx(1);
        elt_compt = distances(1,idx);
        paire = [L_E(elt,1);L_E(elt,2);L_E(elt_compt,1);L_E(elt_compt,2)];
        M = [M,paire];
    end
end

end

