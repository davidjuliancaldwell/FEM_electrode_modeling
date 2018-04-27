function [newPt] = move_electrodes_along_normal(pt,normal,layerFind)

v1 = pt;
v2 = pt+1*normal;

distance = point_to_line_distance(layerFind, v1, v2);

[minimum,idx] = min(distance);

newPt = layerFind(idx,:);

end

