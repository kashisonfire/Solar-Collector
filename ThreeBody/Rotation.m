function Rot = Rotation(A,B)
if A == 1
    Rot = [1 0 0; 0 cos(B) sin(B); 0 -sin(B) cos(B)];
else 
    if A == 2
        Rot = [cos(B) 0 -sin(B); 0 1 0 ; sin(B) 0 cos(B)];
    else
        if A == 3
             Rot = [cos(B) sin(B) 0; -sin(B) cos(B) 0; 0 0 1]
        end
    end
end
end
