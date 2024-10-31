function [center_out, radii_out, center_in, radii_in] = estimate_outer_inner_circle(img,radii,center)
% img: binary image
% center, radii:results of rough circle detection
   
    center = center;
    radii = radii;
    im = img;

    c1 = round(center(1));
    c2 = round(center(2));
    r = round(radii);

    m = 20; % mergin

    if c2-r-m > 0 && c1-r-m > 0 && c2+r+m < size(im,1) && c1+r+m < size(im,2)

        im_window = im(c2-r-m:c2+r+m,c1-r-m:c1+r+m);
        % figure;imagesc(im_window);
        % 
        % figure;imagesc(im_window(:,1:r+m));
        % figure;imagesc(im_window(:,r+m+1:2*r+2*m+1));
        
        %% outer
        im_window_outer = im_window;
        
        %%
        % area 1 (y-axis direction)
        % figure;imagesc(im_window_outer(:,1:r+m));
        for i=1:2*r+2*m+1
            y = im_window_outer(i,1:r+m);
        
            if sum(ischange(double(y),'variance')) >= 3
                id = find(ischange(double(y),'variance') == 1);
                id = id(2);
                im_window_outer(i,id:r+m) = 0;
            end
        end
    
        % figure;imagesc(im_window_outer(:,1:r+m));
    
        % area 2 (y-axis direction)
        % figure;imagesc(im_window_outer(:,r+m+1:2*r+2*m+1));
        for i=1:2*r+2*m+1
            y = im_window_outer(i,r+m+1:2*r+2*m+1);
        
            if sum(ischange(double(y),'variance')) >= 3
                id = find(ischange(double(y),'variance') == 1);
                id = id(sum(ischange(double(y),'variance'))-1)-1;
                im_window_outer(i,r+m+1:r+m+1+id-1) = 0;
            end
        end
    
        % figure;imagesc(im_window_outer(:,r+m+1:2*r+2*m+1));
    
        %%
        % area 3 (x-axis direction)
        % figure;imagesc(im_window_outer(1:r+m,:));
        for i=1:2*r+2*m+1
            x = im_window_outer(1:r+m,i);
        
            if sum(ischange(double(x),'variance')) >= 3
                id = find(ischange(double(x),'variance') == 1);
                id = id(2);
                im_window_outer(id:r+m,i) = 0;
            end
    
        end
    
        % figure;imagesc(im_window_outer(1:r+m,:));
    
        % area 4 (x-axis direction)
        % figure;imagesc(im_window_outer(r+m+1:2*r+2*m+1,:));
        for i=1:2*r+2*m+1
            y = im_window_outer(r+m+1:2*r+2*m+1,i);
        
            if sum(ischange(double(y),'variance')) >= 3
                id = find(ischange(double(y),'variance') == 1);
                id = id(sum(ischange(double(y),'variance'))-1)-1;
                im_window_outer(r+m+1:r+m+1+id-1,i) = 0;
            end
    
        end
    
        % figure;imagesc(im_window_outer(r+m+1:2*r+2*m+1,:));
    
         %% 
        % figure;imagesc(im_window);
    
        im_window_outer = bwareaopen(im_window_outer, 250);
        % figure;imagesc(im_window_outer);
    
        %% inner
        
        %% y-axis direction
        im_window_mask_y = im2bw(zeros(size(im_window)));
        
        % area 1 (y-axis direction)
        % figure;imagesc(im_window(:,1:r+m));
        % figure;imagesc(im_window_mask_y(:,1:r+m));
        for i=1:2*r+2*m+1
            y = im_window(i,1:r+m);
        
            if sum(ischange(double(y),'variance')) >= 2
                id = find(ischange(double(y),'variance') == 1);
                id = id(2);
                im_window_mask_y(i,id:r+m) = 1;
            end
        end
    
        % figure;imagesc(im_window_mask_y(:,1:r+m));
    
        % area 2 (y-axis direction)
        % figure;imagesc(im_window(:,r+m+1:2*r+2*m+1));
        % figure;imagesc(im_window_mask_y(:,r+m+1:2*r+2*m+1));
        for i=1:2*r+21
            y = im_window(i,r+m+1:2*r+2*m+1);
        
            if sum(ischange(double(y),'variance')) >= 2
                id = find(ischange(double(y),'variance') == 1);
                id = id(sum(ischange(double(y),'variance'))-1)-1;
                im_window_mask_y(i,r+m+1:r+m+1+id-1) = 1;
            end
        end
    
        % figure;imagesc(im_window_mask_y(:,r+m+1:2*r+2*m+1));
        
        %% x-axis direction
        im_window_mask_x = im2bw(zeros(size(im_window)));
    
        % area 3 (x-axis direction)
        % figure;imagesc(im_window_mask_x(1:r+m,:));
        for i=1:2*r+2*m+1
            x = im_window(1:r+m,i);
        
            if sum(ischange(double(x),'variance')) >= 2
                id = find(ischange(double(x),'variance') == 1);
                id = id(2);
                im_window_mask_x(id:r+m,i) = 1;
            end
        end
    
        % figure;imagesc(im_window_mask_x(1:r+m,:));
    
        
        % area 4 (x-axis direction)
        % figure;imagesc(im_window_mask_x(r+m+1:2*r+2*m+1,:));
        for i=1:2*r+2*m+1
            y = im_window(r+m+1:2*r+2*m+1,i);
        
            if sum(ischange(double(y),'variance')) >= 2
                id = find(ischange(double(y),'variance') == 1);
                id = id(sum(ischange(double(y),'variance'))-1)-1;
                im_window_mask_x(r+m+1:r+m+1+id-1,i) = 1;
            end
        end
        % figure;imagesc(im_window_mask_x(r+m+1:2*r+2*m+1,:));
    
        %%
        % figure;imagesc(im_window);
        % figure;imagesc(im_window_mask_y);
        % figure;imagesc(im_window.*im_window_mask_y);
        % figure;imagesc(im_window.*im_window_mask_x);
    
        im_window_inner = im_window.*im_window_mask_y.*im_window_mask_x;
        % figure;imagesc(im_window_inner);
    
        im_window_inner = bwareaopen(im_window_inner, 50);
        % figure;imagesc(im_window_inner);
    
        %% outer circle detection
        min_th = round(r*0.8); 
        max_th = r*2;
        x = 0.95; % sensitivity of circular detection, change if necessary
    
        [center_out, radii_out, metric_out] = imfindcircles(im_window_outer,[min_th max_th], Sensitivity=x); 
    
        if length(radii_out) > 1
            p = [];
            for j=1:length(radii_out)
                p_temp = (center_out(j,1) - r)^2 + (center_out(j,2) - r)^2;
                p = [p, p_temp];
            end
        
            [M, I] = min(p);
            
            center_out = center_out(I, 1:2);
            radii_out = radii_out(I)
    
        end
    
        img_circle = draw_circles(rescale(im_window_outer),radii_out,center_out);
        % figure;imagesc(img_circle);
    
        %% inner circle detection
        min_th = round(r*0.8); 
        max_th = r*2;
        x = 0.95; % sensitivity of circular detection, change if necessary
    
        [center_in, radii_in, metric_in] = imfindcircles(im_window_inner,[min_th max_th], Sensitivity=x); 
    
        if length(radii_in) > 1
            p = [];
            for j=1:length(radii_in)
                p_temp = (center_in(j,1) - r)^2 + (center_in(j,2) - r)^2;
                p = [p, p_temp];
            end
        
            [M, I] = min(p);
            
            center_in = center_in(I, 1:2);
            radii_in = radii_in(I)
    
        end
    
        img_circle = draw_circles(rescale(im_window_inner),radii_in,center_in);
        % figure;imagesc(img_circle);
    
        %%
    
        if length(radii_out) + length(radii_in) == 2
            c1_in = center_in(1,1) + (c1-r-m-1);
            c2_in = center_in(1,2) + (c2-r-m-1);
            c1_out = center_out(1,1) + (c1-r-m-1);
            c2_out = center_out(1,2) + (c2-r-m-1);

            c_in = [c1_in, c2_in];
            c_out = [c1_out, c2_out];
        
            center_in = c_in;
            center_out = c_out;
            radii_in = radii_in;
            radii_out = radii_out;
        
        else
            center_in = [];
            center_out = [];
            radii_in = [];
            radii_out = [];
        end
    
        % c_in = [c1_in, c2_in];
        % c_out = [c1_out, c2_out];
        % 
        % center_in = c_in;
        % center_out = c_out;
        % radii_in = radii_in;
        % radii_out = radii_out;
    
    else
        center_in = [];
        center_out = [];
        radii_in = [];
        radii_out = [];

    end