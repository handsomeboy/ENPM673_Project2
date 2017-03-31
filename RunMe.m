%file RunMe.m
%Author: Zejiang
%*@Copyright 2017.2-2017.5 University of Maryland, Zejiang Zeng (zzeng@terpmail.umd.edu)


%%
clear 
close all
% read the image amd remove the noise
OriginalImage=imread('frame38.jpg');
%figure, imshow(OriginalImage);
%remove the noise using median filter
GrayImage=rgb2gray(OriginalImage);
NoNoise=medfilt2(GrayImage);
% figure,imshow(NoNoise);
% Iblur1 = imGaussFilter(NoNoise,2);
% figure,imshow(Iblur1)
%%
% % Image enhencement
% pout_imadjust = imadjust(Iblur1);
% pout_histeq = histeq(Iblur1);
% figure,imshow(pout_imadjust)
% Iblur2 = imGaussFilter(pout_imadjust,2);
% figure,imshow(Iblur2)
%%
% generate horizontal edge emphasis kernel
% h = fspecial('sobel');
% % invert kernel to detect vertical edges
% h = h';
% edge = imfilter(NoNoise,h);
% BW_edge=im2bw(edge,0.7);
% figure,imshow(BW_edge)

%%
BW_edge_canny=edge(NoNoise,'Canny',0.6);
figure,imshow(BW_edge_canny)

%% 
%Mask out the top half of the image
[width,length,~]=size(OriginalImage);
x=[length/6 length*5/6 length*5/6 length/6 length/6];
y=[round(width/2) round(width/2) width width round(width/2)];
Mask=poly2mask(x,y,width,length);
I_Masked = BW_edge_canny.*Mask;
figure,imshow(I_Masked)
%%
% Using Hough transform to find the lines
[H,T,R]=hough(I_Masked);
P  = houghpeaks(H,10,'threshold',ceil(0.1*max(H(:))));
lines = houghlines(I_Masked,T,R,P,'FillGap',20,'MinLength',20);
Pre_lines=lines;
figure,imshow(NoNoise), hold on
max_len = 0;
[~,NumLines]=size(lines);
%%
%Fliter the lines into two grows, right and left lane candidate
for i=NumLines:-1:1
    if abs(lines(i).theta)>70
        lines(i)=[];
    end
end  
[~,NewLines]=size(lines);
p_i=1;
n_i=1;
for j=1:NewLines
    if lines(j).theta>0
        left_lane(p_i)=lines(j);
        p_i=p_i+1;
    else
        right_lane(n_i)=lines(j);
        n_i=n_i+1;
    end
end
%%
% Extrapolate lines
[~,NumLL]=size(left_lane);
for LL_i=1:NumLL
LL_point1(LL_i)=left_lane(LL_i).point1(2);
LL_point2(LL_i)=left_lane(LL_i).point2(2);
end
[~,LL_index_low]=max(LL_point1);
[~,LL_index_high]=min(LL_point2);
LL_point_low=left_lane(LL_index_low).point1;
LL_point_high=left_lane(LL_index_high).point2;
LL_point=[LL_point_low;LL_point_high]';
RL_point1_y=LL_point(2,1);
RL_point2_y=LL_point(2,2);
RL_point1_x=(right_lane(1).rho-RL_point1_y*sind(right_lane(1).theta))/cosd(right_lane(1).theta);
RL_point2_x=(right_lane(1).rho-RL_point2_y*sind(right_lane(1).theta))/cosd(right_lane(1).theta);
RL_point=[RL_point1_x,RL_point2_x;RL_point1_y,RL_point2_y];
%%
for k = 1:NewLines
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

  % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
%plot(xy_long(:,1),xy_long(:,2),'LineWidth',5,'Color','cyan');
plot(LL_point(1,:),LL_point(2,:),'LineWidth',5,'Color','red');
hold on
plot(RL_point(1,:),RL_point(2,:),'LineWidth',5,'Color','red');