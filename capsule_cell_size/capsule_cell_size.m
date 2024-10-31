clear all
close all

%% scale
img_scale = imread('./scale/YG-exp034-aft-vol-2-core-shell-in-media-1-x10_Bottom Slide_D_p00_0_A01f00d4.TIF');
isa(img_scale, 'uint8')
figure;imagesc(img_scale(:,:,1:3));

size(img_scale)

% 445 pixels = 275um
% 1 pixel = 275/445 um
% 1 pixel^2 = (275/445)^2 um^2

unit_conversion = 275/445; % pixel to um

%% set path
path = './input/';
files_d1 = dir([path,'*d1.tif']); % read only GFP images
files_d0 = dir([path,'*d0.tif']); % read only DAPI images
files_d2 = dir([path,'*d2.tif']); % read only RFP images (SYTOX orange)
files_d4 = dir([path,'*d4.tif']); % read only Phase contrast images

savepath ='./output/';

addpath('MyFunctions')

%% prepare output metrics
output_metrics(:,1:5) = 0;

%% kk: Process images one by one
% kk = 1;
for kk = 1:length(files_d1)
% for kk = 4:5

    %% read images
    img_d1 = imread(strcat(path,files_d1(kk).name));
    % max(img_d1, [], "all")
    img_d1 = rescale(img_d1);
    % max(img_d1, [], "all")
    % figure;imagesc(img_d1);
    
    img_d0 = imread(strcat(path,files_d0(kk).name));
    img_d0 = rescale(img_d0);
    % figure;imagesc(img_d0);

    img_d2 = imread(strcat(path,files_d2(kk).name));
    img_d2 = rescale(img_d2);
    % figure;imagesc(img_d2);

    img_d4 = imread(strcat(path,files_d4(kk).name));
    img_d4 = rescale(img_d4);
    % figure;imagesc(img_d4);

    %% GFP image processing: interpolation　of an original image
    sz = size(img_d1);
    xg = 1:sz(1);
    yg = 1:sz(2);

    xq = (0:1/2:sz(1));
    yq = (0:1/2:sz(2));

    [Xg,Yg] = meshgrid(yg,xg);
    [Xq,Yq] = meshgrid(yq,xq);

    img_d1_interp = interp2(Xg, Yg, img_d1, Xq, Yq);
    % figure;imagesc(img_d1_interp);
    size(img_d0)
    size(img_d1_interp)
    
    %% GFP image processing: image preprocessing: denoise by movemean
    img_d1_interp_2 = movmean(movmean(img_d1_interp,20,1),20,2);
    % img_d1_interp_2 = medfilt2(img_d1_interp,[20 20]);
    % figure;imagesc(img_d1_interp_2);

    %% GFP image processing: image preprocessing: calculate gradient using GFP images
    [gx,gy] = gradient(img_d1_interp_2); % 2-dimensional numerical gradient (edge detection)
    grad = sqrt(gx.^2+gy.^2); 
    % figure;imagesc(grad);

    %% GFP image processing: image preprocessing: binalization
    %th = multithresh(grad); % set threshold for binalization
    th = 0.004;
    grad_binary = (grad>th); % create a binary image
    % figure;imagesc(grad_binary);

    %% GFP image processing: rough detection of outer circles: hough transform
    min_th = 70; 
    max_th = 200;
    x = 0.93; % sensitivity of circular detection, change if necessary

    [centers, radiis, metric] = imfindcircles(grad_binary,[min_th max_th], Sensitivity=x); 

    % img_circle = draw_circles(rescale(img_processed),radiis,centers);
    % figure;imagesc(img_circle);

    img_circle = draw_circles(img_d1_interp,radiis,centers);
    % figure;imagesc(img_circle);

    %% GFP image processing: remove large satellites
    % In shells, max intensitiy pixel should be around the edge, whereas it should be around the center as for satellites
    % Remove if the distance between the max intensitiy pixel and the center of the circle is less than 1/2 of the radius 

    [x y] = size(img_d1_interp);
    [xt yt] = meshgrid(1:x, 1:y);

    id = [];

    for i = 1:length(radiis)
        mask_circle = transpose(sqrt((yt-centers(i,1)).^2+(xt-centers(i,2)).^2)<=radiis(i));
        % figure;imagesc(mask_circle);
    
        img_d1_interp_2_masked = img_d1_interp_2.*mask_circle;
        % figure;imagesc(img_d1_interp_2_masked);

        img_d1_interp_2_masked_max = img_d1_interp_2_masked == max(img_d1_interp_2_masked, [], "all");
        % figure;imagesc(img_d1_interp_2_masked_max);

        s = regionprops(img_d1_interp_2_masked_max,'centroid');
    
        if sqrt((s.Centroid(1)-centers(i,1)).^2+(s.Centroid(2)-centers(i,2)).^2) > radiis(i)/2;
            id = [id,i];
        end
    end

    centers = centers(id, 1:2);
    radiis = radiis(id);

    %% show images (with ditected circles)
    img_circle = draw_circles(img_d1_interp,radiis,centers);
    % figure;imagesc(img_circle);

    img_circle_d1 = draw_circles(img_d1,radiis/2,centers/2);
    % figure;imagesc(img_circle);

    img_circle_d0 = draw_circles(img_d0,radiis/2,centers/2);
    % figure;imagesc(img_circle);

    img_circle_d2 = draw_circles(img_d2,radiis/2,centers/2);
    % figure;imagesc(img_circle);

    img_circle_d4 = draw_circles(img_d4,radiis/2,centers/2);
    % figure;imagesc(img_circle);

    %% save images
    
    figure;imagesc(img_circle_d0);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files_d0(kk).name, '_circle', '.pdf'),'ContentType','vector')

    figure;imagesc(img_circle_d1);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files_d1(kk).name, '_circle', '.pdf'),'ContentType','vector')

    figure;imagesc(img_circle_d2);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files_d2(kk).name, '_circle', '.pdf'),'ContentType','vector')

    figure;imagesc(img_circle_d4);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files_d4(kk).name, '_circle', '.pdf'),'ContentType','vector')

    %% Hoechst image processing: interpolation　of an original image
    sz = size(img_d0);
    xg = 1:sz(1);
    yg = 1:sz(2);

    xq = (0:1/2:sz(1));
    yq = (0:1/2:sz(2));

    [Xg,Yg] = meshgrid(yg,xg);
    [Xq,Yq] = meshgrid(yq,xq);

    img_d0_interp = interp2(Xg, Yg, img_d0, Xq, Yq);
    % figure;imagesc(img_d0_interp);
    size(img_d0)
    size(img_d0_interp)

    %% Hoechst image processing: binalization
    img_d0_interp_binary = img_d0_interp > 0.15;
    % figure;imagesc(img_d0_interp_binary);

    % remove small structures
    img_d0_interp_binary_2 = bwareaopen(img_d0_interp_binary,15);
    % figure;imagesc(img_d0_interp_binary_2);

    %% save image
    img_circle_d0_binary_masked = draw_circles(img_d0_interp_binary_2,radiis*0.8,centers);
    % figure;imagesc(img_circle_d0_binary_masked);

    figure;imagesc(img_circle_d0_binary_masked);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files_d0(kk).name, '_binary_masked_circle', '.pdf'),'ContentType','vector')

    %% SYTOX orange image processing: interpolation　of an original image
    sz = size(img_d2);
    xg = 1:sz(1);
    yg = 1:sz(2);

    xq = (0:1/2:sz(1));
    yq = (0:1/2:sz(2));

    [Xg,Yg] = meshgrid(yg,xg);
    [Xq,Yq] = meshgrid(yq,xq);

    img_d2_interp = interp2(Xg, Yg, img_d2, Xq, Yq);
    % figure;imagesc(img_d2_interp);
    size(img_d2)
    size(img_d2_interp)

    %% SYTOX orange image processing: binalization
    img_d2_interp_binary = img_d2_interp > 0.15;
    % figure;imagesc(img_d2_interp_binary);

    % remove small structures
    img_d2_interp_binary_2 = bwareaopen(img_d2_interp_binary,15);
    % figure;imagesc(img_d2_interp_binary_2);

    %% Calculate the areas of Hoechst and SYTOX positive pixels (k: Process shells one by one)
    j = size(output_metrics,1); % Cumulative value of loops

    for k = 1:length(radiis)
    % k = 1;
        %%
        
        mask_circle = transpose(sqrt((yt-centers(k,1)).^2+(xt-centers(k,2)).^2)<=radiis(k)*0.8); % Use 80% of the size of the detected circles to remove edge-covered cells
        % figure;imagesc(mask_circle);

        img_d1_interp_2_masked = img_d1_interp_2.*mask_circle;
        % figure;imagesc(img_d1_interp_2_masked);

        img_d0_interp_binary_2_masked = img_d0_interp_binary_2.*mask_circle;
        % figure;imagesc(img_d0_interp_binary_2_masked);
        % nnz(img_d0_interp_binary_2_masked)*(unit_conversion)^2; % size of Hoechst-positive region in masked images

        img_d2_interp_binary_2_masked = img_d2_interp_binary_2.*mask_circle;
        % figure;imagesc(img_d2_interp_binary_2_masked);
        % nnz(img_d2_interp_binary_2_masked)*(unit_conversion)^2; % size of SYTOX-positive region in masked images

        output_metrics(j+k,1) = kk;
        output_metrics(j+k,2) = k;
        output_metrics(j+k,3) = nnz(img_d0_interp_binary_2_masked)/4*(unit_conversion)^2;
        output_metrics(j+k,4) = nnz(img_d2_interp_binary_2_masked)/4*(unit_conversion)^2;
        output_metrics(j+k,5) = radiis(k)/2*unit_conversion;

        % clearvars output_metrics
    end

    close all
end

%% export table
T = array2table(output_metrics,'VariableNames',{'image_no', 'shell_no', 'cell_area_um2', 'dead_cell_area_um2', 'radii_um'});
T(1,:)=[];

writetable(T, strcat(savepath, "Exp-367-day3-x10.2024-09-09-06-27-26_cell_area_radii_20240920.txt"),'Delimiter','\t')

% histogram(T.cell_area_um2);
% histogram(T.dead_cell_area_um2);

mean(T.radii_um*2)

% ans =
% 
%    69.7071