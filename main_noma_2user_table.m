clear
%close all

ct = 10000;
snr=10000; 

height = 3;
R = 1.9;
snrvec = [100:5 : 130];
M=2;
N=M;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength
neff = 1.4;%low-index materials like Teflon
lambda_g = lambda/neff;
eta = (c/4/pi/f)^2;
D = 5;
D2 = 10;
D1 = 30;

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

        %OMA multiple antenna upper bound
        distance1 = sqrt(loc(1,2).^2+height^2);
        distance2 = sqrt(loc(2,2).^2+height^2);
        oma_r_in(1,i) = 1/2*log2(1+N*M*snr*eta./distance1^2);
        oma_r_in(2,i) = 1/2*log2(1+N*M*snr*eta./distance2^2);

        %analysis for the gap between NOMA and OMA
        diff_in(i) = log2(distance1/distance2)-3;

        %test
        romasum_sim(i) = oma_r_in(1,i) + oma_r_in(2,i);
        romasum(i) = log2(eta*M*snr) -log2(distance1) - log2(distance2)+log2(N);%
        rnomasum_sim(i) = sum(r_all_mul(:,i));
        rnomasum(i) = -log2(distance2^2)+log2(eta*snr)-log2(N);
        diff_ana_test(i) = rnomasum(i) - romasum(i);
        diff_sim_test(i) = sum(r_all_mul(:,i)) - (oma_r_in(1,i) + oma_r_in(2,i));
        %zzz = log2(eta*M*snr) +1/2*log2(N./distance1^2)+1/2*log2(N./distance2^2);
    end    

    %average of average
    r_ave_mul_sum(mi) = mean(sum(r_all_mul,1));% sum rate

    %OMA upper bound
    r_oma(mi) = mean(sum(oma_r_in,1));

    %analysis for the gap of NOMA and OMA
    diff_ave(mi) = mean(diff_in);

end

%plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum,'-o', snrvec,r_ave_mul_sum,'-s', snrvec, r_noma_ana,'-*' )
%plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum,'-o', snrvec,r_ave_mul_sum,'-s')
%plot(snrvec, r_ind,snrvec, r_ind_pin,'-o', snrvec,r_ind_mul,'-s' ) 
%plot(snrvec, r_ave_sum,snrvec, r_ave_pin_sum,'-o',  snrvec,r_oma, '-d', snrvec,r_ave_mul_sum,'-s',snrvec,r_noma_ana,'-s' ) 
%gap_sim = r_ave_mul_sum-r_oma;
%gap_ana = diff_ave;
plot(snrvec,r_ave_mul_sum-r_oma,snrvec,diff_ave)
 