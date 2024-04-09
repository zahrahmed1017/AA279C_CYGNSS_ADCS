function NewInertiaTensor = ParallelAxisTheorem(InertiaTensor, TranslationVector, mass)

    NewInertiaTensor = zeros(3);
    R  = TranslationVector;
    Rnorm = dot(R,R);

    for i = 1:3
        for j = 1:3
            if i == j
                NewInertiaTensor(i,j) = InertiaTensor(i,j) + (mass * (Rnorm - (R(i)*R(j))));
            else
                NewInertiaTensor(i,j) = InertiaTensor(i,j) + (mass * (-R(i) * R(j)));
            end
        end
    end

end