function cygnss = SurfaceDiscretization(cygnss, component)

% Each row of surfaceVertices should be a vertex of the face
% There should be 4 rows for each face (this only works for rectangles).

 faces = fieldnames(cygnss.(component).coord);

 % Iterate through each face of the component
 for i = 1:length(faces)

     faceName = faces{i};
     vertices = cygnss.(component).coord.(faceName);

     T_component_CM = cygnss.(component).origin - cygnss.CM;

     if length(vertices) == 4
         [area, barycenter, normal] = rectangleGeometry(vertices, T_component_CM);

     elseif length(vertices) == 6
         [area, barycenter, normal] = hexagonGeometry(vertices, T_component_CM);
     end

     cygnss.(component).normal.(faceName) = normal;
     cygnss.(component).barycenter.(faceName) = barycenter;
     cygnss.(component).area.(faceName) = area;

end