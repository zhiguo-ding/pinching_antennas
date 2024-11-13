clear
%close all

ct = 50000;
snr=10000; 

height = 3;
R = 1.9;
snrvec = [100:5 : 120];
M=5;
N=M;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength
neff = 1.4;%low-index materials like Teflon
lambda_g = lambda/neff;
eta = (c/4/pi/f)^2;
D = 2;
D2 = 10;
D1 = 20;

%NOMA power allocation 
stepp = 2;
power_temp = flip([1:stepp:stepp*(M-1)+1]);
power = power_temp/sum(power_temp);
for mi = 1 : length(snrvec) 
    snr = 10^(snrvec(mi)/10);

    %conventional antenna
    r_all = zeros(M,ct);
    r_all_pin = zeros(M,ct);
    r_all_mul = zeros(M,ct);
    for i = 1 : ct
        %first general all five users' locations
        Mmax = 5;
        locx = zeros(Mmax,2);
        locx(:,1) = D*rand(Mmax,1)-D/2; %length        
        locx(:,2) = D*rand(Mmax,1)-D/2; %width,  
        locx(Mmax,1) = -D2 + locx(Mmax,1); %best user, shift to the left 
        for m = Mmax-1:-1 : 1
            locx(m,:) = (Mmax-m)*D1 + locx(m,:);%shift to the right and bottom
        end
        loc = [locx(1:M-1,:); locx(Mmax,:)];

        %conventional antenna
        dist = max(1,sqrt(loc(:,1).^2+ loc(:,2).^2 +height^2));
        r_all(:,i) = log2(1+eta*snr./dist.^2)/M;      
        
        %pinching antenna
        dist = max(1,sqrt(loc(:,2).^2 +height^2));
        r_all_pin(:,i) = log2(1+eta*snr./dist.^2)/M;      

        %multiple pinching antennas         
        for mt = 1 : M %for each user, there are M distances to the M pinching antennas
            dist_temp = max(1,sqrt((loc(mt,1)-loc(:,1)).^2 + (loc(mt,2)-0).^2 +height^2));
            theta_temp = 2*pi* (loc(:,1)+100)/lambda_g; %just assume 100 as the location of the feed
            theta_channel = 2*pi*dist_temp/lambda;
            channel_temp(mt) =  abs(sum(1./dist_temp.*exp(-complex(0,1)*(theta_channel+theta_temp))))^2;
        end
            
        for m = 1 : M
            rtemp = zeros(M-m+1,1);
            for mj = m : M
                rtemp(mj-m+1) = log2(1+channel_temp(mj)*power(mj)*eta*snr/M/...
                    (1+channel_temp(mj)*sum(power(mj+1:M))*eta*snr/M));
            end
            r_all_mul(m,i) = min(rtemp);
        end 

    end

    %conventional antenna
    r_ave_sum(mi) = mean(sum(r_all,1));%  sum rate
    for m = 1 : M
        P_all(m,mi) =  length(find(r_all(m,:)<=R))/ct; %individual outage
        r_ind(m,mi) = mean(r_all(m,:)); %individual rate
    end
    %pinching antenna 
    r_ave_pin_sum(mi) = mean(sum(r_all_pin,1)); %sum rate
    for m = 1 : M
        P_all_pin(m,mi) =  length(find(r_all_pin(m,:)<=R))/ct;%individual outage
        r_ind_pin(m,mi) = mean(r_all_pin(m,:)); %individual rate
    end
    

    %average of average
    r_ave_mul_sum(mi) = mean(sum(r_all_mul,1));% sum rate
    for m = 1 : M
        P_all_pin_mul(m,mi) =  length(find(r_all_mul(m,:)<=R))/ct;%individual outage
        r_ind_mul(m,mi) = mean(r_all_mul(m,:)); %individual rate
    end

    %analysis
    r_noma_ana(mi) = -log2(power(2))+log2(D^2/4+height^2+eta*snr*power(2)/2) +2*log2(exp(1))...
        -log2(D^2/4+height^2) - 4/D*log2(exp(1))*height*atan(D/2/height);

end

%plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum,'-o', snrvec,r_ave_mul_sum,'-s', snrvec, r_noma_ana,'-*' )
plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum,'-o', snrvec,r_ave_mul_sum,'-s')
%plot(snrvec, r_ind,snrvec, r_ind_pin,'-o', snrvec,r_ind_mul,'-s' ) 
%plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum,'-o', snrvec,r_ave_mul_sum,'-s',snrvec,r_noma_ana,'-s' ) 
 