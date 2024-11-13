clear
%close all

ct = 500;
snr=100000;
Mvec=[1:5];
D = 20;
height = 3;
R = 1.5;
alp=2;
M=2; 
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2; 
neff = 1.4;%low-index materials like Teflon
lambda_g = lambda/neff;
snrvec = [100:5 : 120];

parfor isnr = 1 : length(snrvec)
    snr=10^(snrvec(isnr)/10);
    R_bound = zeros(ct,1);
    R_mrc = zeros(ct,1);
    R_zf = zeros(ct,1);
    R_pin = zeros(ct,1);
    for ict = 1 : ct

        %users' location
        loc = zeros(M,2);
        loc(1,1) = D*rand(1,1)-D/2; %length
        loc(1,2) = 1*D/6*rand(1,1)+D/3; %width
        loc(2,1) = D*rand(1,1)-D/2; %length
        loc(2,2) = 1*D/6*rand(1,1)-D/2; %width                

        stepsize = lambda/10;% 
        pinvec = [-10*lambda: stepsize : abs(loc(2,1)-loc(1,1))/2]; %no need to search the full segement to avoid duplication
        
        max_value = 0;
        sinr = zeros(length(pinvec),length(pinvec));
        for i = 1 : length(pinvec)
            for j = 1 : length(pinvec)
                %pinching antennas' locations
                pin1 = [loc(1,1)+sign(loc(2,1)-loc(1,1))*pinvec(i) D/3]; %sign is to decide which direction to go
                pin2 = [loc(2,1)-sign(loc(2,1)-loc(1,1))*pinvec(j) -D/3]; 
                %pin1 = [pinvec(i)+loc(1,1) D/3]; 
                %pin2 = [pinvec(j)+loc(2,1) -D/3]; 
        
                %the four channels
                dis11 = sqrt(sum((loc(1,:)-pin1).^2)+height^2);
                dis12 = sqrt(sum((loc(1,:)-pin2).^2)+height^2);
                dis21 = sqrt(sum((loc(2,:)-pin1).^2)+height^2);
                dis22 = sqrt(sum((loc(2,:)-pin2).^2)+height^2);
                %[dis11 dis12 ;dis21 dis22]
                h11 = sqrt(eta)/dis11*exp(-complex(0,1)*2*pi/lambda*dis11)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(1)-pin1(1))); %using pin1(1) as the start of wavguide feed
                h12 = sqrt(eta)/dis12*exp(-complex(0,1)*2*pi/lambda*dis12)*exp(-complex(0,1)*2*pi/lambda_g*(pin2(1)-pin1(1)));
                h21 = sqrt(eta)/dis21*exp(-complex(0,1)*2*pi/lambda*dis21)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(1)-pin1(1)));
                h22 = sqrt(eta)/dis22*exp(-complex(0,1)*2*pi/lambda*dis22)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(2)-pin1(1)));
        
                h1 = [h11 ; h12];
                h2 = [h21; h22];
        
                H = [h1 h2]; 
                %precoding 
                P = inv(H'); %H^H * P = diag
                p1 = P(:,1)/sqrt(P(:,1)'*P(:,1));
                p2 = P(:,2)/sqrt(P(:,2)'*P(:,2));
        
                sinr1 = abs(h1' * p1)^2 * snr;%h11*p1(1) + h12*p1(2);
                sinr2 = abs(h2' * p2)^2 * snr;
                %h21*p2(1) + h22*p2(2);
        
                [sinr(i,j), ind] = min([sinr1 sinr2]);
                if sinr(i,j)>max_value
                    max_value = sinr(i,j);
                    H_max = H;
                end
            end
        end
        
         %performance of multi-waveguide pinching antennas
        R_pin(ict) = log2(1+max(max(sinr)));

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %conventional methods
        % pin1 = [loc(1,1) D/3]; 
        % pin2 = [loc(2,1) -D/3]; 
        pin1 = [0 lambda/4]; 
        pin2 = [0 -lambda/4]; 
        
        %the four channels
        dis11 = sqrt(sum((loc(1,:)-pin1).^2)+height^2);
        dis12 = sqrt(sum((loc(1,:)-pin2).^2)+height^2);
        dis21 = sqrt(sum((loc(2,:)-pin1).^2)+height^2);
        dis22 = sqrt(sum((loc(2,:)-pin2).^2)+height^2);
        %[dis11 dis12 ;dis21 dis22]
        h11 = sqrt(eta)/dis11*exp(-complex(0,1)*2*pi/lambda*dis11)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(1)-pin1(1)));
        h12 = sqrt(eta)/dis12*exp(-complex(0,1)*2*pi/lambda*dis12)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(2)-pin1(1)));
        h21 = sqrt(eta)/dis21*exp(-complex(0,1)*2*pi/lambda*dis21)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(1)-pin1(1)));
        h22 = sqrt(eta)/dis22*exp(-complex(0,1)*2*pi/lambda*dis22)*exp(-complex(0,1)*2*pi/lambda_g*(pin1(2)-pin1(1)));
        
        h1 = [h11 ; h12];
        h2 = [h21; h22];
        
        H = [h1 h2]; 
        P = inv(H'); %H^H * P = diag
        
        %ZF
        p1 = P(:,1)/sqrt(P(:,1)'*P(:,1));
        p2 = P(:,2)/sqrt(P(:,2)'*P(:,2));
        sinr1_zf = abs(h1' * p1)^2 * snr;%h11*p1(1) + h12*p1(2);
        sinr2_zf = abs(h2' * p2)^2 * snr;
        %h21*p2(1) + h22*p2(2);       
        sinr_zf  = min([sinr1_zf sinr2_zf]);
        R_zf(ict) = log2(1+sinr_zf);
        
        %MRC
        sinr1_mrc = h1'*h1/(abs(h1'*h2)^2/(h2'*h2)+1/snr);
        sinr2_mrc = h2'*h2/(abs(h2'*h1)^2/(h1'*h1)+1/snr);
        sinr_mrc  = min([sinr1_mrc sinr2_mrc]);
        R_mrc(ict) = log2(1+sinr_mrc);
        
        %upper bound
        R_bound(ict) = log2(1+ snr*min(abs(h1'*h1),abs(h2'*h2)));
        R1_best = log2(1+ snr*abs(h1'*h1));
        R2_best = log2(1+ snr*abs(h2'*h2));
        % Rmatrix = [rate1 rate2 R_pin; R1_best R2_best R_bound;
        %     log2(1+sinr1_mrc) log2(1+sinr2_mrc) R_mrc;
        %     log2(1+sinr1_zf) log2(1+sinr2_zf) R_zf];
    end
    R_pin_all(isnr) = mean(R_pin)
    R_bound_all(isnr) = mean(R_bound)
    R_mrc_all(isnr) = mean(R_mrc)
    R_zf_all(isnr) = mean(R_zf)
end

plot(snrvec,R_mrc_all,snrvec,R_zf_all,snrvec,R_bound_all,snrvec,R_pin_all)


 %       [R_mrc R_zf R_pin R_bound]
 

