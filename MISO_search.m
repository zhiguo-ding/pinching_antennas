clear
close all

ct = 10000;
snr=100000;
Mvec=[1:5];
D = 10;
height = 3;
R = 1.5;
alp=2;
M=2; 
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2; 

f_cut = 10e9; %cutoff frequency 
lambda_cut = c/f_cut; %cutoff wavelength
lambda_g = lambda/sqrt(1-(lambda/lambda_cut)^2);

snrvec = [100:5 : 120];
snr=10^(snrvec(1)/10);

%users' location
loc = zeros(M,2);
loc(:,1) = D*rand(M,1)-D/2; %length
loc(:,2) = D*rand(M,1)-D/2; %width
        % 
        % loc = zeros(M,2);
        % loc(1,1) = D*rand(1,1)-D/2; %length
        % loc(1,2) = 1*D/6*rand(1,1)+D/3; %width
        % loc(2,1) = D*rand(1,1)-D/2; %length
        % loc(2,2) = 1*D/6*rand(1,1)-D/2; %width   
        
boundary_pin = 10*lambda;
stepsize = lambda/40;
pinvec = [-boundary_pin: stepsize : boundary_pin]; %1st antenna location range, around 1st user

max_value = 0;
sinr = zeros(length(pinvec),length(pinvec));
for i = 1 : length(pinvec)
    for j = 1 : length(pinvec)
        %pinching antennas' locations
        pin1 = [pinvec(i)+loc(1,1) D/3]; 
        pin2 = [pinvec(j)+loc(2,1) -D/3]; 

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

        [sinr(i,j), ind(i,j)] = min([sinr1 sinr2]);
        if sinr(i,j)>max_value
            max_value = sinr(i,j);
            H_max = H;
        end
    end
end

mesh(pinvec,pinvec, sinr)
%performance of multi-waveguide pinching antennas
R_pin = log2(1+max(max(sinr)))
%find the individual data rates
h1x = H_max(:,1); h2x = H_max(:,2); 
Px = inv(H_max'); %H^H * P = diag
p1x = Px(:,1)/sqrt(Px(:,1)'*Px(:,1));
p2x = Px(:,2)/sqrt(Px(:,2)'*Px(:,2));
rate1 = log2(1+snr*abs(h1x' * p1x)^2);
rate2 = log2(1+snr*abs(h2x' * p2x)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%conventional methods
pin1 = [loc(1,1) D/3]; 
pin2 = [loc(2,1) -D/3]; 

%the four channels
dis11 = sqrt(sum((loc(1,:)-pin1).^2)+height^2);
dis12 = sqrt(sum((loc(1,:)-pin2).^2)+height^2);
dis21 = sqrt(sum((loc(2,:)-pin1).^2)+height^2);
dis22 = sqrt(sum((loc(2,:)-pin2).^2)+height^2);
%[dis11 dis12 ;dis21 dis22]
h11 = sqrt(eta)/dis11*exp(-complex(0,1)*2*pi/lambda*dis11);
h12 = sqrt(eta)/dis12*exp(-complex(0,1)*2*pi/lambda*dis12);
h21 = sqrt(eta)/dis21*exp(-complex(0,1)*2*pi/lambda*dis21);
h22 = sqrt(eta)/dis22*exp(-complex(0,1)*2*pi/lambda*dis22);

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
R_zf = log2(1+sinr_zf);

%MRC
sinr1_mrc = h1'*h1/(abs(h1'*h2)^2/(h2'*h2)+1/snr);
sinr2_mrc = h2'*h2/(abs(h2'*h1)^2/(h1'*h1)+1/snr);
sinr_mrc  = min([sinr1_mrc sinr2_mrc]);
R_mrc = log2(1+sinr_mrc);

%upper bound
R_bound = log2(1+ snr*min(abs(h1'*h1),abs(h2'*h2)));
R1_best = log2(1+ snr*abs(h1'*h1));
R2_best = log2(1+ snr*abs(h2'*h2));
Rmatrix = [rate1 rate2 R_pin; R1_best R2_best R_bound;
    log2(1+sinr1_mrc) log2(1+sinr2_mrc) R_mrc;
    log2(1+sinr1_zf) log2(1+sinr2_zf) R_zf]

 %       [R_mrc R_zf R_pin R_bound]
 

