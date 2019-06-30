
function X = conv(x,y)
        X = zeros(1,3);
        if(x<1)
            X(1) = abs(x)+2;
            X(3) = -1;
        elseif(x>400)
            X(1) = 800 - x;
            X(3) = -1;
        else
            X(1) = x;
            X(3) = 1;
        end
        
        %X(1) = mod(-abs(x-1),798);
        %if X(1)>399
         %   X(1) = 798-X(1);
        %end
        %X(1)=X(1)+1;
        X(2) = 1+mod((y-1),10);
        
end
 
 
