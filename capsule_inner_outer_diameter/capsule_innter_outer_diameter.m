clear all
close all

%% scale
img_scale = imread('./scale/YG-exp034-aft-vol-2-core-shell-in-media-1-x10_Bottom Slide_D_p00_0_A01f00d4.TIF');
isa(img_scale, 'uint8')
figure;imagesc(img_scale(:,:,1:3));

% 445 pixels = 275um
unit_conv = 275/445; % (um /1 pixel)

%% set path
path = './input/';

files = dir([path,'*d4.tif']); % read only phase contrast images

savepath ='./output/';

addpath('MyFunctions')

%% prepare output metrics
output_metrics(:,1:9) = 0;

%%
for kk = 1:length(files)
% for kk = 5:6
% kk = 10;
    
    %% read an image
    img = imread(strcat(path,files(kk).name));
    % isa(img, 'uint16')
    % max(img, [], "all")
    img = rescale(img);
    % max(img, [], "all");
    % figure;imagesc(img);

    %% interpolation　of an original image
    sz = size(img);
    xg = 1:sz(1);
    yg = 1:sz(2);

    xq = (0:1/2:sz(1));
    yq = (0:1/2:sz(2));

    [Xg,Yg] = meshgrid(yg,xg);
    [Xq,Yq] = meshgrid(yq,xq);

    img_interp = interp2(Xg, Yg, img, Xq, Yq);
    % figure;imagesc(img_interp);
    % size(img)
    % size(img_interp)

    %% estimate backgroung estimation and interpolation
    ill = imgaussfilt(movmean(movmean(img,200,1),200,2),20);
    % figure;imagesc(ill);
    
    sz = size(img);
    xg = 1:sz(1);
    yg = 1:sz(2);

    xq = (0:1/2:sz(1));
    yq = (0:1/2:sz(2));

    [Xg,Yg] = meshgrid(yg,xg);
    [Xq,Yq] = meshgrid(yq,xq);

    ill_interp = interp2(Xg, Yg, ill, Xq, Yq);
    % figure;imagesc(ill_interp);
   
    %% image preprocessing: normalization of background illumination
    ill_interp = ill_interp / max(ill_interp,[],'all');
    % figure;imagesc(ill_interp);

    img_2 = img_interp./ill_interp;
    % figure;imagesc(img_2);

    %% image preprocessing: denoise by movemean
    % img_3 = medfilt2(img_2,[20 20]);
    % figure;imagesc(img_3);

    img_3 = movmean(movmean(img_2,5,1),5,2);
    % figure;imagesc(img_3);

    %% image preprocessing: calculate gradient
    [gx,gy] = gradient(double(img_3)); % 2-dimensional numerical gradient (edge detection)
    grad = sqrt(gx.^2+gy.^2); 
    % figure;imagesc(grad);
  
    %% image preprocessing: binalization
    %th = multithresh(grad); % set threshold for binalization
    th = 0.008;
    grad_binary = (grad>th); % create a binary image
    % figure;imagesc(grad_binary);

    %% image preprocessing: erosion and dilation
    s =3;
    
    se90 = strel('line',s,90);
    se0 = strel('line',s,0);

    % closing to remove large satellites (dilation->erosion)
    grad_binary_c_1 = imdilate(grad_binary,[se90 se0]);
    % figure;imagesc(grad_binary_c_1);

    grad_binary_c_2 = imerode(grad_binary_c_1,[se90 se0]);
    % figure;imagesc(grad_binary_c_2);

    % Only shell structures remain by opening (erosion->dilation)
    % grad_binary_o_1 = imerode(grad_binary,[se90 se0]);
    grad_binary_o_1 = imerode(grad_binary_c_2,[se90 se0]);
    % figure;imagesc(grad_binary_o_1);

    grad_binary_o_2 = imdilate(grad_binary_o_1,[se90 se0]);
    % figure;imagesc(grad_binary_o_2);

    % remove small structures
    img_processed = bwareaopen(grad_binary_o_2,1500);
    % figure;imagesc(img_processed);

    %% rough detection of outer circles: hough transform
    min_th = 70; 
    max_th = 600;
    x = 0.85; % sensitivity of circular detection, change if necessary

    [centers, radiis, metric] = imfindcircles(img_processed,[min_th max_th], Sensitivity=x); 

    img_circle = draw_circles(rescale(img_processed),radiis,centers);
    % figure;imagesc(img_circle);

    img_circle = draw_circles(img_interp,radiis,centers);
    % figure;imagesc(img_circle);

    %% remove crescent-shaped structures :fill hollow-cores
    % th = multithresh(grad); % set threshold for binalization
    th = 0.009;
    grad_binary_2 = (grad>th); % create a binary image
    % figure;imagesc(grad_binary_2);

    img_filled = imfill(grad_binary_2,'holes');
    % figure;imagesc(img_filled);

    %% remove crescent-shaped structures
    % crescent-shaped structures can be detected as circles, but they are
    % not filled by "imfill"
    id = [];

    for i = 1:length(radiis)
        r=round(centers(i,1));
        c=round(centers(i,2));
    
        if img_filled(c,r) == 1;
            id = [id,i];
        end
    end

        centers = centers(id, 1:2);
    radiis = radiis(id);

    img_circle = draw_circles(rescale(img_processed),radiis,centers);
    % figure;imagesc(img_circle);

    img_circle = draw_circles(img_interp,radiis,centers);
    % figure;imagesc(img_circle);

    %% detection of inner and outer circles
    radiis_in_all = [NaN];
    centers_in_all = [NaN NaN];
    radiis_out_all = [NaN];
    centers_out_all = [NaN NaN];

    % img: binary image
    % center, radii:results of rough circle detection

    % for k = 1:2
    for k = 1:length(radiis)
        [center_out, radii_out, center_in, radii_in] = estimate_outer_inner_circle(img_processed,radiis(k),centers(k,1:2));

        radiis_in_all = cat(1,radiis_in_all,radii_in);
        centers_in_all = cat(1,centers_in_all, center_in);
        radiis_out_all = cat(1,radiis_out_all,radii_out);
        centers_out_all = cat(1,centers_out_all, center_out);
    end

    radiis_in_all = rmmissing(radiis_in_all);
    centers_in_all = rmmissing(centers_in_all);
    radiis_out_all = rmmissing(radiis_out_all);
    centers_out_all = rmmissing(centers_out_all);

    radiis_in_out_all = cat(1,radiis_in_all,radiis_out_all);
    centers_in_out_all = cat(1,centers_in_all,centers_out_all);

    % img_circle = draw_circles(rescale(img_processed),radiis_in_all,centers_in_all);
    % figure;imagesc(img_circle);
    % img_circle = draw_circles(rescale(img_processed),radiis_out_all,centers_out_all);
    % figure;imagesc(img_circle);
    % img_circle = draw_circles(rescale(img_processed),radiis_in_out_all,centers_in_out_all);
    % figure;imagesc(img_circle);

    % img_circle = draw_circles(img_interp,radiis_in_all,centers_in_all);
    % figure;imagesc(img_circle);
    % img_circle = draw_circles(img_interp,radiis_out_all,centers_out_all);
    % figure;imagesc(img_circle);
    % img_circle = draw_circles(img_interp,radiis_in_out_all,centers_in_out_all);
    % figure;imagesc(img_circle);

    img_circle_1 = draw_circles(img,radiis_in_all/2,centers_in_all/2);
    % figure;imagesc(img_circle_1);
    img_circle_2 = draw_circles(img,radiis_out_all/2,centers_out_all/2);
    % figure;imagesc(img_circle_2);
    img_circle_3 = draw_circles(img,radiis_in_out_all/2,centers_in_out_all/2);
    % figure;imagesc(img_circle_3);

    %% save images
    % imwrite(uint8(img_circle*255), strcat(savepath, files(kk).name,'_circle','.tiff'))
    
    figure;imagesc(img_circle_1);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files(kk).name, '_circle_in', '.pdf'),'ContentType','vector')

    figure;imagesc(img_circle_2);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files(kk).name, '_circle_out', '.pdf'),'ContentType','vector')


    figure;imagesc(img_circle_3);
    colorbar
    axis image
    f = gcf;
    exportgraphics(gcf,strcat(savepath, files(kk).name, '_circle_in_out', '.pdf'),'ContentType','vector')

    %% save metrics
    temp(:,1) = repelem(kk,length(radiis_out_all));
    temp(:,2) = [radiis_out_all/2];
    temp(:,3) = [radiis_out_all/2*unit_conv*2];
    temp(:,4:5) = [centers_out_all/2];
    temp(:,6) = [radiis_in_all/2];
    temp(:,7) = [radiis_in_all/2*unit_conv*2];
    temp(:,8:9) = [centers_in_all/2];
    output_metrics = cat(1,output_metrics, temp);

    clearvars temp

end

%% export table

T = array2table(output_metrics,'VariableNames',{'image_no','outer_radius_pixel',...
    'outer_diameter_um', 'center_outer_1','center_outer_2','inner_radius_pixel',...
    'inner_diameter_um', 'center_inner_1','center_inner_2'});
T(1,:)=[];

c = split(pwd,"/");
MyDir = char(c(length(c)));
writetable(T, strcat(savepath, MyDir, '.txt'),'Delimiter','\t')

% figure;histogram(T.outer_diameter_um);
% f = gcf;
% exportgraphics(gcf,strcat(savepath, MyDir, '_outer_histgram', '.pdf'),'ContentType','vector')
% 
% figure;histogram(T.inner_diameter_um);
% f = gcf;
% exportgraphics(gcf,strcat(savepath, MyDir, '_inner_histgram', '.pdf'),'ContentType','vector')
% 
% std(T.outer_diameter_um)/mean(T.outer_diameter_um)*100
% std(T.inner_diameter_um)/mean(T.inner_diameter_um)*100
