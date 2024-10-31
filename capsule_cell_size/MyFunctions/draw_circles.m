function [output] = draw_circles(img,radii,centers) % radii, centers: array (output of imfindcircles)
    %% draw circles
    img_circle = img;

    for k = 1:length(radii)     
        theta=0:1:360;
        r=round(centers(k,1) + radii(k)*sin(theta));
        c=round(centers(k,2) + radii(k)*cos(theta));
        
        r(r < 1) = 1;
        c(c < 1) = 1;
        r(r > size(img,2)) = size(img,2);
        c(c > size(img,1)) = size(img,1);
    
        for j=1:length(r)
            img_circle(c(j),r(j))=2;
        end
    end

    output = img_circle;

end