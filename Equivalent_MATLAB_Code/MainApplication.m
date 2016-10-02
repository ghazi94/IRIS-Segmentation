%function to detect the pupil  boundary
%it searches a certain subset of the image
%with a given radius range(rmin,rmax)
%around a 10*10 neighbourhood of the point x,y given as input
%INPUTS:
%im:image to be processed
%rmin:minimum radius
%rmax:maximum radius
%x:x-coordinate of centre point
%y:y-coordinate of centre point
%OUTPUT:
%cp:centre coordinates of pupil(x,y) followed by radius
%Author:Anirudh S.K.
%Department of Computer Science and Engineering
%Indian Institute of Techology,Madras
function [ci]=search(im,rmin,rmax,x1,y1,option)
rows=size(im,1);
cols=size(im,2);
sigma=0.5;%(standard deviation of Gaussian)
R=rmin:rmax;
maxrad=zeros(rows,cols);
maxb=zeros(rows,cols);
for i=(x1-5):(x1+5)
for j=(y1-5):(y1+5)
        [b,r,blur]=partiald(im,[i,j],rmin,rmax,0.5,600,option);
        maxrad(i,j)=r;
        maxb(i,j)=b;
    end
end
B=max(max(maxb));
[X,Y]=find(maxb==B);
% disp('checkpoint1_search');
r1=maxrad(X,Y);
ci=[X,Y,r1];



        