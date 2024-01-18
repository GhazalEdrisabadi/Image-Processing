computer = imread('SchoolOfComputer.jpg');
landscape = imread('Landscape.jpg');
me = imread('Me.jpg');
figure(1);
imshow(computer);
figure(2);
imshow(landscape);
figure(3);
imshow(me);

% sobel algorithm to detect edges
clc; clear; close all; warning off;
c = imread('Me.jpg');
c = im2double(c);
c = imnoise(c, 'gaussian', 0, 0.001);

subplot(2,2,1);
imshow(c);
title('orginal image');

[height width R] = size(c);

for i = 2 : height-1
    for j = 2 : width-1
        Dx(i,j) = [c(i+1,j-1) - c(i-1,j-1)] + 2*[c(i+1,j) - c(i-1,j)] + [c(i+1,j+1) - c(i-1,j+1)];
        Dy(i,j) = [c(i-1,j+1) - c(i-1,j-1)] + 2*[c(i,j+1) - c(i,j-1)] + [c(i+1,j+1) - c(i+1,j-1)];
        S(i,j) = sqrt(Dx(i,j)^2 + Dy(i,j)^2);
        if Dx(i,j) < 1
            Dx(i,j) = 0;
        else
            Dx(i,j) = 1;
        end
        if Dy(i,j) < 1
            Dy(i,j) = 0;
        else
            Dy(i,j) = 1;
        end
    end
end

subplot(2,2,2);
imshow(Dx,[]);
title('horizontal edges detected');

subplot(2,2,3);
imshow(Dy,[]);
title('vertical edges detected');

for i = 1 : 255
    for j = 1 : 255
       if (S(i,j) < 1)
            S(i,j) = 0;
       else
           S(i,j) = 1;
       end
    end
end
subplot(2,2,4);
imshow(S,[]);
title('final image');

% roberts algorithm to detect edges
clc; clear; close all; warning off;
I = imread('Me.jpg');
I = imnoise(I, 'gaussian', 0, 0.001);
I = im2double(I);

subplot(2,2,1);
imshow(I);

[height width R] = size(I);
for i = 2 : height-1
    for j = 2 : width-1
        R(i,j) = abs(I(i+1,j+1)-I(i,j)) + abs(I(i+1,j)-I(i,j+1));
        Z(i,j) = abs(I(i+1,j+1)-I(i,j));
        X(i,j) = abs(I(i+1,j)-I(i,j+1));
    end
end

for i = 1 : height-1
    for j = 1 : width-1
        if (R(i,j) < 0.25)
            R(i,j) = 0;
        else
            R(i,j) = 1;
        end
        if (Z(i,j) < 0.25)
            Z(i,j) = 0;
        else
            Z(i,j) = 1;
        end
        if (X(i,j) < 0.25)
            X(i,j) = 0;
        else
            X(i,j) = 1;
        end
    end
end

subplot(2,2,2);
imshow(Z,[]);

subplot(2,2,3);
imshow(X,[]);

subplot(2,2,4);
imshow(R,[]);

%------------ matlab functions to edge detecting ------------

I = rgb2gray(imread("Me.jpg"));
subplot(1,2,1);
imshow(I);  
title("Gray Scale Image");  

% Sobel Edge Detection  
S = edge(I, 'Sobel');  
subplot(1,2,2);
imshow(S);  
title("Sobel");

% Robert Edge Detection  
R = edge(I, 'Roberts');
subplot(1,2,2); 
imshow(R);  
title("Robert");

% Canny Edge Detection  
C = edge(I, 'Canny');
subplot(1,2,2); 
imshow(C);  
title("Canny");

% Zerocross Edge Detection  
Z = edge(I, 'zerocross');
subplot(1,2,2);
imshow(Z);  
title("Zerocross");

%------------ matlab functions to make noisy images ------------

% read image from directory
I = imread('Me.jpg');
%subplot(1,2,1);
imshow(I);
title("original image");

% add gaussian noise => avr=0.1 , var=0.02
G = imnoise(I, 'gaussian', 0.1, 0.02);
%subplot(1,2,2);
figure(3);
imshow(G);
title("gaussian noise");

% add salt & pepper noise => avr=0.1 , var=0.02
SP = imnoise(I,'salt & pepper',0.02);
%subplot(1,2,2);
figure(4);
imshow(SP);
title("salt & pepper noise");

% add poisson noise
P = imnoise(I, 'poisson');
subplot(1,2,2);
imshow(P);
title("poisson noise");

% add speckle noise => var=0.02
SPC = imnoise(I, 'speckle', 0.02);
subplot(1,2,2);
imshow(SPC);
title("speckle noise");

%------------ matlab functions to remove noises ------------

% remove gaussian noise : first we make noisy gray image, then remove
% noises with wiener2 filtering
I = rgb2gray(imread("Me.jpg"));
J = imnoise(I,'gaussian',0,0.025);
subplot(1,2,1);
imshow(J);
title("noisy image");
K = wiener2(J,[5 5]);
subplot(1,2,2);
imshow(K);
title("image without noise");
