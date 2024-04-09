function [area, barycenter, normal] = rectangleGeometry(vertices)

    % Calculate the vectors representing two adjacent sides of the rectangle
    side1 = vertices(2,:) - vertices(1,:);
    side2 = vertices(4,:) - vertices(1,:);
    
    % Calculate the lengths of these sides
    lengthSide1 = norm(side1);
    lengthSide2 = norm(side2);
    
    % Area of the rectangle is the product of the lengths of its sides
    area = lengthSide1 * lengthSide2;
    
    % Calculate barycenter by averaging the coordinates
    barycenter = mean(vertices, 1);
    
    % Calculate the normal vector by taking the cross product of two sides
    normal = cross(side1, side2);
    % Normalize the normal vector to have unit length
    normal = normal / norm(normal);
    

end