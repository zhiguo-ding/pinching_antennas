clear
%close all

ct = 10000;
snr=10000;
Mvec=[1:5];
D = 40;
height = 3; %d 
snrvec = [100:5 : 120];
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;

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

        %disk boundary, upper bound for conventional antenna
        mtemp=0;
        loc_up = zeros(M,2);
        while mtemp<M
            loc_up_temp = D*rand(1,2)-D/2;
            if sqrt(loc_up_temp*loc_up_temp')<D/2
                mtemp = mtemp + 1;
                loc_up(mtemp,:) = loc_up_temp;
            end
        end
        %upper bound
        dist_up = max(1,sqrt(loc_up(:,1).^2+ loc_up(:,2).^2 +height^2));
        r_all_up(:,i) = log2(1+eta*snr./dist_up.^2)/M; 

    end

    %conventional antenna
    r_ave_sum(mi) = mean(sum(r_all,1));%  sum rate 
    %conventional antenna upper bound
    r_ave_sum_up(mi) = mean(sum(r_all_up,1));%  sum rate 
    %pinching antenna 
    r_ave_pin_sum(mi) = mean(sum(r_all_pin,1)); %sum rate 

    %analytical results
    a=height^2; b=eta*snr;
    tau2 = D/2*log2(D^2/4+a+b);
    term1 = tau2 - log2(exp(1))*D + 2*log2(exp(1))*sqrt(a+b)*atan(D/2/sqrt(a+b));
    a=height^2; b=0;
    tau2 = D/2*log2(D^2/4+a+b);
    term2 = tau2 - log2(exp(1))*D + 2*log2(exp(1))*sqrt(a+b)*atan(D/2/sqrt(a+b));
    r_ana_pin(mi) = 2*(term1 -term2)/D;

    %analytical results - upper bound
    a=height^2+eta*snr;
    temp1 = D^2/4*log(D^2/4+a)-D^2/4+a*log((D^2/4+a)/a);
    a=height^2;
    temp2 = D^2/4*log(D^2/4+a)-D^2/4+a*log((D^2/4+a)/a);
    r_all_up_ana(mi) = 4/D^2*log2(exp(1))*(temp1-temp2);

    %high SNR approximation
    r_ana_pin_app(mi) =  log2(D^2/4+height^2+eta*snr)+2*log2(exp(1)) ...
        - log2(D^2/4+height^2)- 4/D*log2(exp(1))*height*atan(D/2/height);
    r_all_up_app(mi) = log2(D^2/4+height^2+eta*snr) + log2(exp(1)) ...
        -log2(D^2/4+height^2) -4/D^2*height^2*log2((D^2/4+height^2)/height^2);

    %test
    diff1(mi) = r_ana_pin_app(mi)-r_all_up_app(mi);
    x=D/height/2;
    diff2(mi) = log2(exp(1))-2./x*log2(exp(1)).*atan(x) +1./x.^2.*log2(1+x.^2);

end
 
plot(snrvec, r_ave_sum,snrvec, r_ave_sum_up,'-x',snrvec,r_all_up_ana,'-s',snrvec, r_all_up_app, '-d',snrvec, r_ave_pin_sum,'-o', snrvec, r_ana_pin, snrvec, r_ana_pin_app, '-*' ) 
 