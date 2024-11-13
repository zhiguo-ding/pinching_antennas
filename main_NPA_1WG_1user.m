%one waveguide, N pinching antennas, one user
clear
%close all

ct = 5000;
Mvec=[1:5];
D = 10;
height = 3; %d 
snrvec = [100:5 : 120];
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
neff = 1.4;%low-index materials like Teflon
lambda_g = lambda/neff;
eta = (c/4/pi/f)^2;
theta_eps = 0.001;
N=3;%number of pinching antennas
guide_loc = lambda/2;

for mi = 1 : length(snrvec) 
    snr = 10^(snrvec(mi)/10);
    %NOMA power allocation 
    power = [4/5 1/5];
    %conventional antenna
    r_all = zeros(M,ct);
    r_all_pin = zeros(M,ct);
    r_all_mul = zeros(M,ct);
    for i = 1 : ct
        loc = zeros(M,2);
        loc(:,1) = D*rand(M,1)-D/2; %length        
        loc(:,2) = D*rand(M,1)-D/2; %width,         


        %conventional antenna
        dist = max(1,sqrt(loc(:,1).^2+ loc(:,2).^2 +height^2));
        r_all(:,i) = log2(1+eta*snr./dist.^2)/M;      
        
        %pinching antenna
        dist = max(1,sqrt(loc(:,2).^2 +height^2));
        r_all_pin(:,i) = log2(1+eta*snr./dist.^2)/M;   

        %multiple pinching antennas
        %find the location of the first one
        stepsize = lambda_g/100;
        find_if =0;
        loc_new = [loc(1,1) 0];% P_1^pin is closest to the first user
        while find_if == 0
            loc_new = [loc_new(1,1)+stepsize 0]; %move it to the right
            theta_n = 2*pi*(loc_new(1,1)+D/2)/lambda_g; %phase shift in the waveguide
            theta_channel = 2*pi/lambda*sqrt(sum((loc(1,:)-loc_new).^2)+height^2);%phase shift of the channel
            delta_theta = (theta_n+theta_channel)/2/pi;
            diff_theta = abs(round(delta_theta)-delta_theta);
            if diff_theta<theta_eps
                find_if = 1; % got the first antennas
            end
        end
        loc_pin(1) = loc_new(1,1);

        % find the n-th location 
        for n = 2 : N
            find_if =0;
            loc_new = [loc_pin(n-1)+guide_loc 0];% P_1^pin is closest to the first user
            while find_if == 0
                loc_new = [loc_new(1,1)+stepsize 0]; %move it to the right
                theta_n = 2*pi*(loc_new(1,1)+D/2)/lambda_g; %phase shift in the waveguide
                theta_channel = 2*pi/lambda*sqrt(sum((loc(1,:)-loc_new).^2)+height^2);%phase shift of the channel
                delta_theta = (theta_n+theta_channel)/2/pi;
                diff_theta = abs(round(delta_theta)-delta_theta);
                if diff_theta<theta_eps
                    find_if = 1; % got the first antennas
                end
            end
            loc_pin(n) = loc_new(1,1);
        end
        %find the data rate
        loc_pin_2D = [loc_pin' zeros(N,1)]; %NX2
        distn = sqrt(sum((loc(1,:)-loc_pin_2D).^2,2)+height^2); %from nth antenna to the first user
        thetan = 2*pi*(loc_pin_2D(:,1)+D/2)/lambda_g; %phase shift in the waveguide; theta_n
        channel_m = eta*abs(sum( exp(-1i*2*pi/lambda*distn)./distn.*exp(-1i*thetan) ) )^2;        
        r_mul_pin(i) = log2(1+channel_m*snr/N);

        %upper bound
        dist_shortest = sqrt(sum((loc(1,:)-[loc(1,1) 0]).^2,2)+height^2);  
        r_mul_pin_upp(i) = log2(1+N*snr*eta/dist_shortest^2);
    end

    %conventional antenna
    r_ave_sum(mi) = mean(sum(r_all,1));%  sum rate  
    %pinching antenna 
    r_ave_pin_sum(mi) = mean(sum(r_all_pin,1)); %sum rate 
    %user multiple antennas
    r_mul_pin_sum(mi) = mean(r_mul_pin);
    %upper bound
    r_mul_pin_sum_up(mi) = mean(r_mul_pin_upp);
 


end
 
plot(snrvec, r_ave_sum  , snrvec, r_ave_pin_sum,snrvec,r_mul_pin_sum,snrvec,r_mul_pin_sum_up) 
%plot(snrvec,r_mul_pin_sum,snrvec,r_mul_pin_sum_up) 
 