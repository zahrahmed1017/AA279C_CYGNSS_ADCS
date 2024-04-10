function [area, barycenter, normal] = hexagonGeometry(vertices, T_wrtCM)
    % Ensure the vertices are in order (either clockwise or counter-clockwise)
    % vertices is a 6x3 matrix where each row represents a vertex
    
    % Calculate the area using the Shoelace formula
    y = vertices(:,2);
    z = vertices(:,3);
    area = 0.5 * abs(sum(y(1:end-1).*z(2:end)) + y(end)*z(1) - sum(y(2:end).*z(1:end-1)) - y(1)*z(end));
    
    % Calculate barycenter by averaging the coordinates
    barycenter = mean(vertices, 1) + T_wrtCM;
    
    % Calculate two vectors from the hexagon's edges
    vector1 = vertices(2,:) - vertices(1,:);
    vector2 = vertices(4,:) - vertices(1,:);
    
    % Calculate the normal vector by taking the cross product of two edge vectors
    normal = cross(vector1, vector2);
    % Normalize the normal vector to have unit length
    normal = normal / norm(normal);
end