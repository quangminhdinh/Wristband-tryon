clear all; clc; close all;                                 %clear stored variables, command windows, close any open windows
objects = imaqfind;                                        %find video input objects in memory
delete(objects);

%config camera stream
vid=videoinput('winvideo',1);
set(vid,'ReturnedColorspace','rgb')
triggerconfig(vid, 'manual');
start(vid);

while true
    figure(1);
    for idx=1:6
        subplot(2, 3, idx);
        I=imread(string(idx) + ".png");
        M = repmat(all(~I,3),[1 1 3]); %mask black parts
        I(M) = 255; %turn them white
        imshow(I);
        title(idx);
    end
    sgtitle('Please pick a wristband'); 
    
    inp = input("Enter a number from 1-6, 0 to exit: ");
    while inp < 0 || inp > 6
        inp = input("Please enter a number from 1-6: ");
    end
    if inp == 0
    break
    end
    IM1=getsnapshot(vid);                                      %get snapshot from the webcam video and store to IM1 variable
    
    orgwb = imread(string(inp) + ".png");
    wb=orgwb;
    fgh = figure(2);
    while true
        if ~ishghandle(fgh)
            break
        end
        IM2=getsnapshot(vid);                                  %get snapshot of test image and store to variable IM2
        
        
        IM3 = IM1 - IM2;                                                            %subtract Backround from Image
        IM3 = rgb2gray(IM3);                                                        %Converts RGB to Gray
        lvl = graythresh(IM3);                                                      %find the threshold value using Otsu's method for black and white
        
        IM3 = im2bw(IM3, lvl);                                                      %Converts image to BW, pixels with value higher than threshold value is changed to 1, lower changed to 0
        IM3 = bwareaopen(IM3, 10000);
        IM3 = imfill(IM3,'holes');
        IM3 = imerode(IM3,strel('disk',15));                                        %erode image
        IM3 = imdilate(IM3,strel('disk',20));                                       %dilate iamge
        IM3 = medfilt2(IM3, [5 5]);                                                 %median filtering
        IM3 = bwareaopen(IM3, 10000);                                               %finds objects, noise or regions with pixel area lower than 10,000 and removes them
        
        REG=regionprops(IM3,'all');                                                 %calculate the properties of regions for objects found 
        CEN = cat(1, REG.Centroid);                                                 %calculate Centroid
        [B, L, N, A] = bwboundaries(IM3,'noholes');                                 %returns the number of objects (N), adjacency matrix A, object boundaries B, nonnegative integers of contiguous regions L
        
        for k =1:length(B)                                                      %for the given object k            
            BND = B{k};                                                         %boundary set for object
            BNDx = BND(:,2);                                                    %Boundary x coord
            BNDy = BND(:,1);                                                    %Boundary y coord
                
            if (length(B) == 1)
                pkoffset = CEN(:,2);                                              %Calculate peak offset point from centroid
                [pks,locs] = findpeaks(-BNDy,'minpeakheight',-pkoffset);         %find peaks in the boundary in y axis with a minimum height greater than the peak offset
                n_X = BNDx(locs);
            end  
        
        end
        if length(CEN) > 1
            try
            centroids = round(CEN);
    %        start = centroids(1, :);
            if exist('pks','var') && length(pks) > 3
                dvec = [n_X(3), -pks(3)] - centroids(1, :);
                start = dvec * -0.5 + centroids(1, :);
                
                distance = norm(dvec);
                wb = imresize(orgwb, distance / 264);
    
                %CosTheta = max(min(dot(dvec, [0, -1])/norm(dvec),1),-1);
                %angle = real(acosd(CosTheta));
                % a = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
                angle = atan2d(-dvec(1), -dvec(2));
                wb = imrotate(wb,angle,'bilinear','crop');
            end
            if ~exist('start','var')
                start = centroids(1, :);
            end
            ystart = start(2) - fix(size(wb, 1) / 2);
            xstart = start(1) - fix(size(wb, 2) / 2);
    
            rmask = find(wb(:, :, 1));
            gmask = find(wb(:, :, 2));
            bmask = find(wb(:, :, 3));
            wbmask = union(rmask, gmask);
            wbmask = union(wbmask, bmask);
            xmask = fix(wbmask/size(wb, 1));
            ymask = rem(wbmask, size(wb, 1));
    
            
            catch
            end
        end
        if exist('ystart','var')
            try
            %IM2(ystart + ymask, xstart + xmask,:) = wb(ymask, xmask+1,:);
            for cc=0:2
                cmask = cc * size(IM2, 1) * size(IM2, 2) + (xstart + xmask) * size(IM2, 1) + ystart + ymask;
                nwbmask = cc * size(wb, 1) * size(wb, 2) + wbmask;
                IM2(cmask) = wb(nwbmask);
            end
            catch
            end
        end
        figure(2);imshow(IM2);
    end
end
delete(vid);