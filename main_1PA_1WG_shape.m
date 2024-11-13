clear
%close all

ct = 10000;
snr=10000;
Mvec=[1:5];
D = 10;
height = 3; %d 
snrvec = [100:5 : 120];
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;
D_leng =4*D;

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
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length        
        loc(:,2) = D*rand(M,1)-D/2; %width,      

        %conventional antenna
        dist = max(1,sqrt(loc(:,1).^2+ loc(:,2).^2 +height^2));
        r_all(:,i) = log2(1+eta*snr./dist.^2)/M;      
        
        %pinching antenna
        dist = max(1,sqrt(loc(:,2).^2 +height^2));
        r_all_pin(:,i) = log2(1+eta*snr./dist.^2)/M;   
 

    end

    %conventional antenna
    r_ave_sum(mi) = mean(sum(r_all,1));%  sum rate  
    %pinching antenna 
    r_ave_pin_sum(mi) = mean(sum(r_all_pin,1)); %sum rate 
 
end
 
%plot(snrvec, r_ave_sum) 
plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum ) 
 