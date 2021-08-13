clc; clear; close;
%handling only single channel images at present
filename = "2.jpeg";
%Size of the psf
psf_size = 25;
%Type of the psf
psf_type = "gaussian";
%use this if the psf needs another parameter
psf_second_param = 5;

Img_blur = imread(filename);
%from here your image will be a double
Img_blur = double(Img_blur);

imshow(uint8(Img_blur));
title('Our blurred Image');
pause(2);
%remove any dotted noise
noise_range = 10;
Img_blur = removeDottednoise(Img_blur, noise_range);
%End of initial perp
%%
%lets make make a psf estimate, which is close to the original
[~,~,d] = size(Img_blur);
img_ch = Img_blur;
im_in = Img_blur;
%%
%define all the parameters below
%for try2
lambda=0;
n=30;%no. of iterations play a key role (as the 'n' increases deblurr increases but there is a tradeoff between ringing and deblurr)

%for try6
S = numel(Img_blur);
alpha=0.7;

%for try7
t = 20;

choice = "n";
while (choice ~= "y")
    %weight of the image
    edge_weight = input("Enter a suitable value for edge weight (0-10)(can be float) : ");
    imgEdge = edge(im2double(double(rgb2gray(uint8(Img_blur)))),'sobel',edge_weight);
    se = strel('disk',3);
    imgEdge = bwareaopen(imgEdge,3);
    M = 1-double(imdilate(imgEdge,se));
    imshow(M);
    title('Weight');
    drawnow;
    pause(2);
    disp("Note: The black area should be less and the white area more.");
    disp("     Also it should have minimal black dots and have a clear edge.");
    disp("Are you satisfied with the edge Map (y/n)");
    choice = input("","s");
end
imgEdge = bwareaopen(imgEdge,8);
M1 = 1-double(imdilate(imgEdge,se));
M = M1;
%%
for e = 1:d
    Img_blur = img_ch(:,:,e);
    
    if psf_type == "disk"
        PSF = fspecial(psf_type,psf_size);
    else
        PSF = fspecial(psf_type,psf_size,psf_second_param);
    end

    %try2
    inp=Img_blur;
    J1=Img_blur;
    J2=zeros(size(J1));
    P1=PSF;
    P2=zeros(size(P1));
    sizeI=size(J1(:));
    sizeI=sizeI(1);
    sizeP=size(P1(:));
    sizeP=sizeP(1);
    g=zeros(sizeI,2);
    p=zeros(sizeP,2);
    idx = repmat({':'},[1 length(size(J1))]);
    jn=5;
    bn=5;

    Y=max(J1+lambda*(J1-J2),0);
    B=max(P1+lambda*(P1-P2),0);
    sumPsf=sum(B(:));
    B=B/(sumPsf + (sumPsf==0)*eps);
    for k=(lambda+1):(lambda+n)

        for j_it=1:jn
            
            CC=richardlucy(Y,B,inp,idx);
            temp = isnan(CC);
            if sum(temp(:)) > 0
                break
            end
            J2=J1;
            H=psf2otf(P1,size(J1));
            J1=max(Y.*(real(ifftn(conj(H).*CC))),0);
            g=[J1(:)-Y(:) g(:,1)]; 

            if(k>2)
                lambda=((g(:,1).'*g(:,2))/(g(:,2).'*g(:,2)));  %(2)
                lambda=max(min(lambda,1),0);                   %(3)
            end
            Y=max(J1+lambda*(J1-J2),0); %step 1 in algorithm 
%             imshow(uint8(J1));
%             drawnow;
        end
        for b_it=1:bn
            
            CC=richardlucy(Y,B,inp,idx);
            temp = isnan(CC);
            if sum(temp(:)) > 0
                break
            end
            P2=P1;
            H=fftn(J2);
            H( isnan(H) ) = temp( isnan(H) );
            P1=max(B.*(otf2psf(conj(H).*CC,size(P1))),0);
            sumPsf=sum(P1(:));
            P1=P1/(sumPsf + (sumPsf==0)*eps) ;    %normalising psf so that summation of pixel values of psf == 1
            p=[P1(:)-B(:) p(:,1)];
            if(k>2)                                             
                lambda=((p(:,1).'*p(:,2))/(p(:,2).'*p(:,2)));  %same repeats for psf estimation
                lambda=max(min(lambda,1),0);
            end
            B=max(P1+lambda*(P1-P2),0);
            sumPsf=sum(B(:));
            B=B/(sumPsf + (sumPsf==0)*eps);
%             imshow(P1,[],'InitialMagnification','fit');
%             drawnow;
        end
        J1 = removeDottednoise(J1,10);
    end

    img1 = uint8(J1);
    PSF = P1;
    imshow(img1);
    title('Result from try2');
    pause(1);
    %%
    %try6
    %setup
    P=psf2otf(PSF,size(Img_blur));
    %start
    Y=fftn(Img_blur);
    %fn (||AP-Y||.^2)./2 for each pixel
    %make an intial guess for A and other algorithm paramters
    A=Y;
    %since we used dot operator the following operation is done for every pixel
    %the above cost fn is there to just to have a record of the current cost
    total_cost=10000;
    old_cost=20000;

    while (old_cost-total_cost)>0.001
        old_cost = total_cost;
        A=max(A-alpha.*(A.*P-Y).*P,0);
        current_cost=(((real(ifftn(A.*P))-real(ifftn(Y))).^2)./2);
        total_cost=abs(sum(current_cost(:)));
        total_cost=total_cost/S;
    %     disp(total_cost);
    %     imshow(uint8(real(ifftn(A))));
    %     drawnow;
    end
    %convert to spatial domain
    A=real(ifftn(A));
    A = removeDottednoise(uint8(A),30);
    img2 = uint8(A);
    imshow(img2);
    title('Result from try6');
    pause(1);
    %%
    %try7
    img = im2double(Img_blur);
    [delx, dely] = gradient(img);
    %initially
    imgdeblurred = img;
    %the algorithm
    for i = 1:t
        %edge map
        %for x axis
        x = max(conj(delx).*((delx.*imgdeblurred) - (delx.*img)),0);
        %for y axis
        y = max(conj(dely).*((dely.*imgdeblurred) - (dely.*img)),0);
        %edge map
        Epl = 2.*(x+y).*M;

        %for the function
        %1st part
        one_part = imgdeblurred./(1+lambda.*Epl);
        %2nd part
        two_part = imfilter(img./imfilter(imgdeblurred,PSF,'conv'), PSF);
        imgdeblurred = one_part.*two_part;
    end
    imgdeblurred = uint8(imgdeblurred);
    img3 = uint8(removeDottednoise(imgdeblurred, 20));
    imshow(img3);
    title('Result from try7');
    pause(2);
    %%
    %our result
    img_ch(:,:,e) = uint8(0.15*img1 + 0.85*img2 + 0.05*img3);
end
%%
imwrite(uint8(img_ch),"deblurred.jpeg");
subplot(1,2,1);
imshow(uint8(im_in));
title('Our blurred Image');
subplot(1,2,2);
imshow(uint8(img_ch));
title('Our recovered Image');
set(gcf, 'Position', get(0, 'Screensize'));
load gong.mat;
sound(y);
%%
%functions we used
function imgFixed = removeDottednoise(img, clearing_level)
    % img : the input image
    % clearing_level : specifies the noise pixel level
    % imgFixed : the output img
    if (~(exist('clearing_level','var')))
        clearing_level = 0;
    end
    [~, ~, d] = size(img);
    
    channel = img;
    fixed = channel;

    for i = 1:d
        fixed(:, :, i) = medfilt2(channel(:, :, i), [3 3]);
    end

    imgFixed = channel;
    %change this range to remove more noise
    for i = 0:clearing_level
        for j = 1:d
        noiseFree = imgFixed(:, :, j);
        noiseFreevalue = fixed(:, :, j);
        noiseImage = (channel(:, :, j) == i | channel(:, :, j) == (255-i) );
        noiseFree(noiseImage) = noiseFreevalue(noiseImage);
        imgFixed(:, :, j) = noiseFree;
        end
    end
end

function [imgratio] = richardlucy(Y,B,inp,idx)
    H=psf2otf(B,size(Y));
    Y=fftn(Y);
    reblurr=real(ifftn(H.*Y));
    otp=double(inp)./reblurr;
    otp=otp(idx{:});
    imgratio=fftn(otp);
end
