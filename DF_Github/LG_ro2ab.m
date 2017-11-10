function [a,b] = LG_ro2ab(rho)
    b = asin(rho(3));
    
    x = rho(1);
    y = rho(2);
    
    if(x>0.001)
        if(y>0)%第一象限 + +
            a = atan((rho(2))/rho(1));
        else%第四象限 + -
            a = atan((rho(2))/rho(1));
        end
    elseif(x<-0.001)
        if(y>0)%第二象限 -, +
            a = pi - atan(abs(rho(2)/rho(1)));
        else%第三象限 - -
            a = pi + atan(abs(rho(2)/rho(1)));
        end
    else
        a = pi/2;
    end
end