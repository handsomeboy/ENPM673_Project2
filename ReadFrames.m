%*@File ReadFrames.m
%*@Author Zejiang Zeng
%*@Copyright 2017.2-2017.5 University of Maryland, Zejiang Zeng (zzeng@terpmail.umd.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read image from video
clear
v=VideoReader('../Data/project_video.mp4');
for n=490:500
   filename=strcat('frame',num2str(n),'.jpg');
   b=read(v,n);
   imwrite(b,filename);
end
%%