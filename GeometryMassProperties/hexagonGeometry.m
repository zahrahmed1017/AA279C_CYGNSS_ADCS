function [area, barycenter, normal] = hexagonGeometry(vertices)
    % Ensure the vertices are in order (either clockwise or counter-clockwise)
    % vertices is a 6x3 matrix where each row represents a vertex
    
    % Calculate the area using the Shoelace formula
    x = vertices(:,1);
    y = vertices(:,2);
    area = 0.5 * abs(sum(x(1:end-1).*y(2:end)) + x(end)*y(1) - sum(x(2:end).*y(1:end-1)) - x(1)*y(end));
    
    % Calculate barycenter by averaging the coordinates
    barycenter = mean(vertices, 1);
    
    % Calculate two vectors from the hexagon's edges
    vector1 = vertices(2,:) - vertices(1,:);
    vector2 = vertices(4,:) - vertices(1,:);
    
    % Calculate the normal vector by taking the cross product of two edge vectors
    normal = cross(vector1, vector2);
    % Normalize the normal vector to have unit length
    normal = normal / norm(normal);
end