function quad = Targ2Quadrants(startPos, targetPos)

target_xy = targetPos - startPos;


for i = 1:size(target_xy,1)
    if (sign(target_xy(i,1)) == sign(target_xy(i,2)))
        if (target_xy(i,1) > 0) quad(i) = 1;
        else quad(i) = 3;
        end
    else if(target_xy(i,1) > 0) quad(i) = 4;
        else quad(i) = 2;
        end
    end
end