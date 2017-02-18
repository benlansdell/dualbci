function oct = Targ2Octants(startPos, targetPos)

	target_xy = targetPos - startPos;
	for i = 1:size(target_xy,1)
	    if (sign(target_xy(i,1)) == sign(target_xy(i,2)))
	        if (target_xy(i,1) > 0) 
	        	quad(i) = 1;
	        	if (target_xy(i,2)< target_xy(i,1))
		        	oct(i) = 1;
				else
					oct(i) = 2;
				end
	        else
	        	quad(i) = 3;
	        	if (target_xy(i,2)< target_xy(i,1))
		        	oct(i) = 6;
				else
					oct(i) = 5;
				end
	        end
	    else if(target_xy(i,1) > 0) 
	    		quad(i) = 4;
	        	if abs(target_xy(i,2)) < abs(target_xy(i,1))
		        	oct(i) = 8;
				else
					oct(i) = 7;
				end
	        else 
	        	quad(i) = 2;
	        	if abs(target_xy(i,2)) < abs(target_xy(i,1))
		        	oct(i) = 4;
				else
					oct(i) = 3;
				end
	        end
	    end
	end
end