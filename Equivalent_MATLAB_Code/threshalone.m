tic
cd(fileparts(mfilename('fullpath'))) %Sets matlab directory to local directory
I=imread('image1.jpg');
I=imcomplement(imfill(imcomplement(I),'holes'));
cd(fileparts(mfilename('fullpath'))) %Sets matlab directory to local directory
I=imresize(I,0.5);
figure, imshow(uint8(I));
csvwrite('1.csv',I);
I=csvread('1.csv');
rmin=16;
rmax=30;
scale=1;
rows=size(I,1);
cols=size(I,2);
[X,Y]=find(I<110);
%keyboard();
%Generates a column vector of the image elements
%that have been selected by tresholding;one for x coordinate and one for y
s=size(X,1);
for k=1:s %
    if (X(k)>rmin)&(Y(k)>rmin)&(X(k)<=(rows-rmin))&(Y(k)<=(cols-rmin))
            A=I((X(k)-1):(X(k)+1),(Y(k)-1):(Y(k)+1));
            M=min(min(A));
            %this process scans the neighbourhood of the selected pixel
            %to check if it is a local minimum
           if I(X(k),Y(k))~=M
              X(k)=NaN;
              Y(k)=NaN;
           end
    end
end

v=find(isnan(X));

%%keyboard();
X(v)=[];
Y(v)=[];
hold on;
%keyboard();
%deletes all pixels that are NOT local minima(that have been set to NaN)
index=find((X<=rmin)|(Y<=rmin)|(X>(rows-rmin))|(Y>(cols-rmin)));
X(index)=[];
Y(index)=[];
%keyboard();
%This process deletes all pixels that are so close to the border 
%that they could not possibly be the centre coordinates.
%plot(Y, X,'g+', 'MarkerSize', 3);
N=size(X,1);
%recompute the size after deleting unnecessary elements
maxb=zeros(rows,cols);
maxrad=zeros(rows,cols);
%keyboard();
%defines two arrays maxb and maxrad to store the maximum value of blur
%for each of the selected centre points and the corresponding radius
for j=1:N
    [b,r,blur]=partiald(I,[X(j),Y(j)],rmin,rmax,'inf',600,'pupil');%coarse search
    maxb(X(j),Y(j))=b;
    maxrad(X(j),Y(j))=r;
end
[x1,y1]=find(maxb==max(max(maxb)));
%radi=maxrad(x,y);
%r1=maxrad(x1,y1);
%ci=[x1,y1,r1];
%ci=[x1,y1,radi];
%plot(y1, x1,'r+', 'MarkerSize', 6);
cp=search(I,rmin,rmax,x1,y1,'pupil');%fine search
%finds the maximum value of blur by scanning all the centre coordinates
cp=cp/scale
%plot(cp(2), cp(1),'b+', 'MarkerSize', 15);
%the function search searches for the centre of the pupil and its radius
%by scanning a 10*10 window around the iris centre for establishing 
%the pupil's centre and hence its radius
ci=search(I,round(1.9*cp(3)),round(5*cp(3)),cp(1)*scale,cp(2)*scale,'iris')
%Ref:Daugman's paper that sets biological limits 
%on the relative sizes of the iris and pupil
out=drawcircle(I,[cp(1),cp(2)],cp(3),600);
out2=drawcircle(out,[ci(1),ci(2)],ci(3),600);
imshow(uint8(out2));
%out=drawcircle(out,[cp(1),cp(2)],cp(3),600);
hold off;
toc