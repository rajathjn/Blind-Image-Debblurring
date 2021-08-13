clc; clear; close;

%filename = input("Enter the image name : ",'s');
filename = "1_deblurred.jpg";
Image = imread(filename);

blur = fspecial('gaussian',5,1);

Img_blur = real(ifftn( fftn(Image).*psf2otf(blur,size(Image)) ));

imwrite(uint8(Img_blur),"blurred.jpg");