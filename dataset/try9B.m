clc; clear; close;

imgs = 1:54;
%Add your image names here
filename = string(imgs)+'.jpg';
%Make a range of size estimate
%Higher range higher chance you'll find your psf
%But the compute type will also increase
sizes = 2:30;
results = 'result' + string(imgs);
for res = results
    mkdir(res);
end

%index of images you want to run on
for i = 50:54
    log = fopen('info'+string(i)+'.txt','w+');
    fprintf(log,'psfSize,cost,outVarlap1,outVarlapnorm1,outVarlap2,outVarlapnorm2,outMean,outMeannorm\n');   
    disp("testing initialised..");
    tic
    %handling only single channel images at present
    Img_blur = double(rgb2gray(imread(filename(i))));
    disp(filename(i));
%     imshow(uint8(Img_blur));
%     title('Our blurred Image');
%     pause(2);

    %remove any dotted noise
    noise_range = 10;
    Img_blur = removeDottednoise(Img_blur, noise_range);
    %End of initial perp
    %%
    psf_type = "disk";
    %use this if the psf needs another parameter
    psf_second_param = 0;
    %%
    for c = sizes
        if psf_type == "disk"
            PSF = fspecial(psf_type,c);
        else
            PSF = fspecial(psf_type,c,psf_second_param);
        end
        %%
        %define all the parameters below
        %for try2
        lambda=0;
        n=30;%no. of iterations play a key role (as the 'n' increases deblurr increases but there is a tradeoff between ringing and deblurr)

        %for try6
        S = numel(Img_blur);
        alpha=1;
        laplacian1 = [0 1 0 ; 1 -4 1 ; 0 1 0];
        laplacian2 = [0 -1 0; -1 4 -1; 0 -1 0];

        %for try7
        t = 200;
        %%
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
        B=B/(sumPsf + (sumPsf==0)*eps) ;
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
        %         imshow(uint8(J1));
        %         drawnow;
            end
            for b_it=1:bn
                CC=richardlucy(Y,B,inp,idx);
                temp = isnan(CC);
                if sum(temp(:)) > 0
                    break
                end
                P2=P1;
                H=fftn(J2);
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
        %         imshow(P1,[],'InitialMagnification','fit');
        %         drawnow;
            end
            J1 = removeDottednoise(J1,10);
        end

        img1 = uint8(J1);
        PSF = P1;
        % imshow(img1);
        % title('Result from try2');
        % pause(1);
        %%
        %try6
        %setup
        P=psf2otf(PSF,size(Img_blur));
        %start
        Y=fftn(Img_blur);
        Yprime = real(ifftn(Y));
        %fn (||AP-Y||.^2)./2 for each pixel
        %make an intial guess for A and other algorithm paramters
        A=Y;
        total_cost=10000;
        old_cost=20000;

        while (old_cost-total_cost)>0.001
            old_cost = total_cost;
            A = max(A-alpha.*(A.*P-Y).*P,0);
            current_cost = (((real(ifftn(A.*P))-Yprime).^2)./2);
            total_cost = abs(sum(current_cost(:)));
            total_cost = total_cost/S;
        end
        %convert to spatial domain
        A=real(ifftn(A));
        A = removeDottednoise(uint8(A),30);
        img2 = uint8(A);

        psfSize = c;

        temp = imfilter(double(img2),laplacian1,'conv');
        outVarlap1 = sum( temp(:) );
        outVarlapnorm1 = outVarlap1/S;
        temp = imfilter(double(img2),laplacian2,'conv');
        outVarlap2 = sum( temp(:) );
        outVarlapnorm2 = outVarlap2/S;

        temp = double(img2);
        outMean = mean(temp(:));
        outMeannorm = outMean/S;

        fprintf(log,string(psfSize)+','+string(total_cost)+','+string(outVarlap1)+','+string(outVarlapnorm1)+','+string(outVarlap2)+','+string(outVarlapnorm2)+','+string(outMean)+','+string(outMeannorm)+'\n');

        %display results
        % imshow(img2);
        % title('Result from try6');
        % pause(1);
        %our result
        img = 0.4*img1 + 0.6*img2;
        imwrite(img,'result'+string(i)+'\size_'+string(c)+'.jpg');

    end
    fclose(log);
    toc
end

disp("Done");
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
