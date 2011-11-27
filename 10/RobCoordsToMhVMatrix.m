function MhV = RobCoordsToMhVMatrix(P)

    MhV = ...
        [ anglesToMtx(P(4:6,1)), P(1:3,1); ...
          zeros(1,3),            1 ];

end