%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SIMULATOR 5GPencil  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the simulator synthesizes the beams, aimed at the deployment spots, to serve 
% users based on the knowledge of their position and a certain value of 
% localization uncertainty. Then, after the deployment, it computes the 
% SINR and throughput over the users and the contribution in terms of 
% EMF over the serving cell, i.e., the one surrounded by 6 cells.

% Authors: Simone Rossetti, Sara Saida 
% Date: 21.01.2021
% Version: 1.1
%% Parameters setup

close all
clear all
clear
clf

tic

%System Parameters
%%%Reference txPowerDBm: https://www.itu.int/en/ITU-T/Workshops-and-Seminars/20171205/Documents/S3_Christer_Tornevik.pdf
%%%Reference Am, SLAv: https://ieeexplore.ieee.org/document/7465794
%%%Refence GdB: https://halberdbastion.com/products/antenna-catalogue/hamilton-35-ghz-8x8-mimo-panel-antenna use-case Range: 3600 to 3800 MHz
%%%Reference Gtx: https://ieeexplore.ieee.org/document/7997079

UE_h = 1.5; %user height (m)
BS_h = 15; %BS height (m)
effective_h = BS_h-UE_h; %effective height (m)
txPowerDBm = 53; % Total transmit power in dBm for each array 
txPower = (10.^((txPowerDBm-30)/10)); % Convert dBm to W
Am = 25; % Front-back ratio (dB) 
SLAv = 20; %Side-lobe level limit
GdB = 15;  %dBi Gain of the Small Cell
G = 10^(GdB/10); %from dBi to linear Gain
Gtx=3; %antenna element Tx Gain (dBi)
fq=3.7 ; %Carrier Frequency 3.7GHz
Z0=377; %Free space wave impedence
angle_min=3; %minimum value of the angle in degrees

% Define receiver parameters using Table 8-2 (b) of Report ITU-R
% M.[IMT-2020.EVAL] https://www.itu.int/md/R15-SG05-C-0057/es
bw = 80e6; % 80 MHz bandwidth
rxNoiseFigure = 5; % dB  
rxNoisePowerdB = -174 + 10*log10(bw) + rxNoiseFigure;
rxNoisePower = 10^(rxNoisePowerdB/10);

% Define array size (64 elements)
nrow = 8;
ncol = 8;
txPowerSE = txPower/(nrow*ncol); %tx power from each antenna element
EIRP = txPowerSE*G; % Equivalent Isotropic Radiated Power
EIRP_db = 10*log10(EIRP); %Equivalent Isotropic Radiated Power in dB
PG_4pi= (EIRP)/(4*pi); 
Gbf=10*log10(nrow*ncol); %maximum beamforming Gain dB


%% Loading variables from layout generator

load('layout.mat'); %load layout (BS, cell, sectors)
load('deploy.mat'); %load deployment spot
load('measure.mat'); %load measure spot
load('useful_spots.mat'); %load all the spots inside the cells
measure_on=0; %1 to enable computation of e.m field, 0 to turn off
fixed_emf = 0; %1 to enable computation of emf when no beamforming is used
plot_on=0; %1 to enable plot, 0 to disable
NLOS_on=0; %1 for NLOS, 0 for LOS
int_intrasec=0; %1 to consider intrasector interference, 0 to disable
fig_style=0; %to enable plot of user tolerance areas

%% User Generation

accuracy=2; %diameter (m) of area around deployment spot in which the UE is deployed
users{1,size(deploy_spot,2)}=[];

for i=1:size(deploy_spot,2)
    for k=1:length(deploy_spot{1,i})
        pos=[deploy_spot{1,i}(k,1), deploy_spot{1,i}(k,2)]; %position of deployment spot
        c=0;
        while c<1 
            ru= (accuracy/2)*sqrt(rand(1,1)); %to have uniform distribution
            theta=2*pi*rand(1,1);
            xr=pos(1,1)+ru*cos(theta);
            yr=pos(1,2)+ru*sin(theta);
            pr=[xr,yr];
            if isinterior(sector{1,i}, pr)==1 %check if it stands inside sector
                pos_user=pr;
                c=1;
                if isinterior(secside_r{1,i},pos_user)==1 %right side
                    users{1,i} = [users{1,i}; pos_user,1];
                else %left side
                    users{1,i} = [users{1,i}; pos_user,0];
                end
            end
        end
    end
end

%% Compute Steer, Tilt, Path loss for Deploy spot and Users

deploy_spot_comp{size(deploy_spot,1),size(deploy_spot,2)}=[]; %to store all deployment spot data
users_comp{size(users,1),size(users,2)}=[]; %to store all users data
for j=1:size(deploy_spot,1)
    for i=1:size(deploy_spot,2)
        for k=1:length(deploy_spot{j,i})
                       
            pos=[deploy_spot{j,i}(k,1), deploy_spot{j,i}(k,2)]; %position of current deploy spot
            center=[BS(j,1), BS(j,2)]; %coordinates of current BS
            dist_2D=norm(pos-center); %2D distance between deploy spot and BS
            dist_eff_ds=sqrt((effective_h)^2 + (dist_2D)^2); %effective distance
            PL_ds = 21*log10(dist_eff_ds) + 32.4 + 20*log10(fq); %Path loss LOS dB
            if NLOS_on==1 %apply NLOS model
                PL_nlos_ds=35.3*log10(dist_eff_ds)+22.4+21.3*log10(fq)-0.3*(UE_h-1.5); %Path loss NLOS dB
                PL_ds=max(PL_nlos_ds,PL_ds);
            end
            %compute tilt angle for deploy spot
            tilt=atan2d(effective_h,dist_2D);
            %compute steer angle for deploy spot
            steer=evalsteer(i,pos,center);
            deploy_spot_comp{j,i}=[deploy_spot_comp{j,i}; deploy_spot{j,i}(k,:), tilt, steer, PL_ds, dist_eff_ds]; %store all info about ds
            
            if j==1 %computation for users only in serving cell
                pos_user=[users{j,i}(k,1), users{j,i}(k,2)]; %position of user
                dist_2D=norm(pos_user-center); %2D distance between user and BS
                dist_eff_user=sqrt((effective_h)^2 + (dist_2D)^2); %effective distance
                PL = 21*log10(dist_eff_user) + 32.4 + 20*log10(fq); %Path loss LOS dB
                if NLOS_on==1 %apply NLOS model
                    PL_nlos=35.3*log10(dist_eff_user)+22.4+21.3*log10(fq)-0.3*(UE_h-1.5); %Path loss NLOS dB
                    PL=max(PL_nlos,PL);
                end
                %compute tilt angle for user
                tilt_u=atan2d(effective_h,dist_2D);
                %compute steer angle for user
                steer_u=evalsteer(i,pos_user,center);
                users_comp{j,i}=[users_comp{j,i}; users{j,i}(k,:), tilt_u, steer_u, PL, dist_eff_user]; %store all info about users
            end
        end
    end
end

save('deploy_comp.mat','deploy_spot_comp'); %all deployment spots with their infos
save('users_comp.mat', 'users_comp'); %all users with their infos
%% Beam Synthesize 
%%%Reference az3dB=65°, el3dB=65°: https://ieeexplore.ieee.org/document/8904367
%%%Reference az3dB=90°, el3dB=90°: https://halberdbastion.com/products/antenna-catalogue/hamilton-35-ghz-8x8-mimo-panel-antenna , use-case Range: 3600 to 3800 MHz

angle_3dB{length(cell),size(sector,2)}=[];
tp_ds{1,size(deploy_spot_comp,2)}=[];
tp_user{1,size(users_comp,2)}=[];
delta_tp{1,size(users_comp,2)}=[];
delta_dist{1,size(users_comp,2)}=[];
dist_user{1,size(users_comp,2)}=[];

for j=1:size(deploy_spot_comp,1)
    for i=1:size(deploy_spot_comp,2)
        for k=1:length(deploy_spot_comp{j,i})
            center=[BS(j,1), BS(j,2)]; %coordinates of current BS
            pos=[deploy_spot{j,i}(k,1), deploy_spot{j,i}(k,2)]; %position of current deploy spot
            dis=norm(pos-center); %2D distance between worst spot and BS
            rc=accuracy/2; %radius of circular area 

            %compute az3dB in order to cover the area of accuracy in which
            [xout,yout]=circcirc(center(1,1),center(1,2),dis,pos(1,1),pos(1,2),rc); 
            p1=[xout(1,1),yout(1,1)];
            p2=[xout(1,2),yout(1,2)];
            dp1=norm(center-p1); dp2=norm(center-p2); %distance points-center    
            dpp=norm(p1-p2);
            dp1_eff=sqrt((dp1^2)+(effective_h^2)); %effective distance from center
            dp2_eff=sqrt((dp2^2)+(effective_h^2));
            %cos(alpha)=((b^2)+(c^2)-10^2)/(2*b*c)) with a,b,c side of triangle
            az3dB=acosd(((dp1_eff^2)+(dp2_eff^2)-(dpp^2))/(2*dp1_eff*dp2_eff));
            if az3dB < angle_min
                az3dB = angle_min;
            end

            %%az3dB=65; %Reference at the beginning of this section 
            %%az3dB=90; %Reference at the beginning of this section 

            %%%compute el3dB
            theta = 0 : 0.01 : 2*pi;
            xc = pos(1,1) + rc * cos(theta);
            yc = pos(1,2) + rc * sin(theta);
            %compute nearest and farthest point of the circumference from BS
            near=dis;
            far=near;
            for n=1:length(xc)
                point=[xc(1,n) yc(1,n)];
                tmp=norm(center-point);
                if tmp<near
                    near=tmp;
                    nearestp=point;
                end
                if tmp>far
                    far=tmp;
                    farthest=point;
                end
            end
            dpp=norm(farthest-nearestp); %point-point distance
            near_eff=sqrt((near^2)+(effective_h^2)); %effective distance from center
            far_eff=sqrt((far^2)+(effective_h^2));
            %cos(alpha)=((b^2)+(c^2)-10^2)/(2*b*c)) with a,b,c side of triangle
            el3dB=acosd(((near_eff^2)+(far_eff^2)-(dpp^2))/(2*near_eff*far_eff));
            if el3dB < angle_min
                el3dB = angle_min;
            end

            %%el3dB=65; %Reference at the beginning of this section 
            %%el3dB=7; %Reference at the beginning of this section 
            
            angle_3dB{j,i}=[angle_3dB{j,i}; az3dB, el3dB]; %to store all angles of each deploy_spot
        end
    end
end

%% compute SINR and tp over deploy spot

for i=1:size(deploy_spot_comp,2)
    j=1;
    for k=1:length(deploy_spot_comp{j,i})
        I_ds_co_so=[];
        I_ds_co_s=[];
        I_ds_c_s=[];
        pos=[deploy_spot_comp{j,i}(k,1),deploy_spot_comp{j,i}(k,2)]; %position of deployment spot
        
        %Define antenna radiation pattern
        A_H = evalAH(deploy_spot_comp, j, i, k, deploy_spot_comp{j,i}(k,5), angle_3dB{j,i}(k,1), Am ); %azimuth radiation pattern
        A_H=10*log10(A_H); %from linear to dB
        A_V = evalAV(deploy_spot_comp, j, i, k, deploy_spot_comp{j,i}(k,4), angle_3dB{j,i}(k,2), SLAv ); %elevation radiation pattern
        A_V=10*log10(A_V); %from linear to dB
        A_tx=A_H+A_V; %overall antenna radiation pattern in dB
        A_tx(A_tx<-Am) = -Am; %overall radiation pattern 
        %compute Beamforming Gain
        BF = 10*log10(real(sinc(((deploy_spot_comp{j,i}(k,5)-deploy_spot_comp{j,i}(k,5))/(1.13*angle_3dB{j,i}(k,1))).^2))) +...
            10*log10(real(sinc(((deploy_spot_comp{j,i}(k,4)-deploy_spot_comp{j,i}(k,4))/(1.13*angle_3dB{j,i}(k,2))).^2))) + Gbf;
        BF=real(BF);
        PL_ds=deploy_spot_comp{j,i}(k,6); %path loss in dB
        P_ds_co = txPowerDBm - PL_ds + A_tx + Gtx + BF; %compute desired signal power 
        P_ds_co = 10^(P_ds_co/10); %from dB to linear
        
        %compute intra sector interference
        for m=1:length(deploy_spot_comp{j,i})
            if m~=k 
                A_H = evalAH(deploy_spot_comp, j, i, m, deploy_spot_comp{j,i}(k,5), angle_3dB{j,i}(m,1), Am ); %azimuth radiation pattern
                A_H=10*log10(A_H); %from linear to dB
                A_V = evalAV(deploy_spot_comp, j, i, m, deploy_spot_comp{j,i}(k,4), angle_3dB{j,i}(m,2), SLAv ); %elevation radiation pattern
                A_V=10*log10(A_V); %from linear to dB
                A_tx=A_H+A_V; %overall antenna radiation pattern in dB
                A_tx(A_tx<-Am) = -Am; %overall radiation pattern
                %compute Beamforming Gain
                BF = 10*log10(sinc(((deploy_spot_comp{j,i}(k,5)-deploy_spot_comp{j,i}(m,5))/(1.13*angle_3dB{j,i}(m,1))).^2)) +...
                    10*log10(sinc(((deploy_spot_comp{j,i}(k,4)-deploy_spot_comp{j,i}(m,4))/(1.13*angle_3dB{j,i}(m,2))).^2)) + Gbf;
                BF=real(BF);
                %compute interference
                I = txPowerDBm - PL_ds + A_tx + Gtx + BF; %intra sector interference
                I = 10^(I/10); %from dB to linear
                I_ds_co_so = [I_ds_co_so, I];
            end
        end
        
        %compute inter sector interference (same cell)
        for l=1:size(deploy_spot_comp,2)
            if l~=i 
                for m=1:length(deploy_spot_comp{j,l})
                    center=[BS(j,1),BS(j,2)]; %coordinates of current BS
                    steer = evalsteer(l, pos, center); %computes steer
                    A_H = evalAH(deploy_spot_comp, j, l, m, steer, angle_3dB{j,l}(m,1), Am ); %azimuth radiation pattern
                    A_H=10*log10(A_H); %from linear to dB
                    A_V = evalAV(deploy_spot_comp, j, l, m, deploy_spot_comp{j,i}(k,4), angle_3dB{j,l}(m,2), SLAv ); %elevation radiation pattern
                    A_V=10*log10(A_V); %from linear to dB
                    A_tx=A_H+A_V; %overall radiation pattern in dB
                    A_tx(A_tx<-Am) = -Am; %overall radiation pattern
                    %compute Beamforming Gain
                    BF = 10*log10(real(sinc(((steer-deploy_spot_comp{j,l}(m,5))/(1.13*angle_3dB{j,l}(m,1))).^2))) +...
                        10*log10(real(sinc(((deploy_spot_comp{j,i}(k,4)-deploy_spot_comp{j,l}(m,4))/(1.13*angle_3dB{j,l}(m,2))).^2))) + Gbf;
                    BF=real(BF);
                    %compute interference
                    I = txPowerDBm - PL_ds + A_tx + Gtx + BF; %inter sector interference
                    I = 10^(I/10); %from dB to linear
                    I_ds_co_s = [I_ds_co_s, I];
                end
            end
        end
        
        %compute inter cell interference
        for n=2:size(deploy_spot_comp,1)
            for l=1:size(deploy_spot_comp,2)
                for m=1:length(deploy_spot_comp{n,l})
                    center=[BS(n,1),BS(n,2)]; %coordinates of current BS
                    d=norm(pos-center); %2D distance between deployment spot and BS
                    d_eff=sqrt((effective_h)^2 + (d)^2); %effective distance
                    tilt=atan2d(effective_h,d); %computes tilt angle
                    steer = evalsteer(l, pos, center); %computes steer angle
                    A_H = evalAH(deploy_spot_comp, n, l, m, steer, angle_3dB{n,l}(m,1), Am ); %azimuth radiation pattern
                    A_H=10*log10(A_H); %from linear to dB
                    A_V = evalAV(deploy_spot_comp, n, l, m, tilt, angle_3dB{n,l}(m,2), SLAv ); %elevation radiation pattern
                    A_V=10*log10(A_V); %from linear to dB
                    A_tx=A_H+A_V; %overall radiation pattern in dB
                    A_tx(A_tx<-Am) = -Am; %overall radiation pattern
                    %compute Beamforming Gain
                    BF = 10*log10(real(sinc(((steer-deploy_spot_comp{n,l}(m,5))/(1.13*angle_3dB{n,l}(m,1))).^2))) +...
                        10*log10(real(sinc(((tilt-deploy_spot_comp{n,l}(m,4))/(1.13*angle_3dB{n,l}(m,2))).^2))) + Gbf;
                    BF=real(BF);
                    %compute path loss
                    PL_ext=21*log10(d_eff) + 32.4 + 20*log10(fq);
                    if NLOS_on==1
                        PL_nlos_ext=35.3*log10(d_eff)+22.4+21.3*log10(fq)-0.3*(UE_h-1.5);
                        PL=max(PL_nlos_ext,PL_ext);
                    end
                    %compute interference
                    I = txPowerDBm - PL_ext + A_tx + Gtx + BF; %inter sector interference
                    I = 10^(I/10); %from dB to linear
                    I_ds_c_s = [I_ds_c_s, I];
                end
            end
        end
        
        %compute SINR 
        if int_intrasec==1
            SINR = P_ds_co / (sum(I_ds_co_so)+sum(I_ds_co_s)+sum(I_ds_c_s)+rxNoisePower);
        else
            SINR = P_ds_co / (sum(I_ds_co_s)+sum(I_ds_c_s)+rxNoisePower);
        end
        instant_ach_rate = real(log2(1+(SINR)));
        tp_ds{j,i}(k)=bw*instant_ach_rate; %throughput computed over deployment spot
        
    end
end


%% compute SINR and tp over users

for i=1:size(users_comp,2)
    j=1;
    for k=1:length(users_comp{j,i})
        I_u_co_so=[];
        I_u_co_s=[];
        I_u_c_s=[];
        pos=[users_comp{j,i}(k,1),users_comp{j,i}(k,2)]; %position of user
        
        %Define antenna radiation pattern
        A_H = evalAH(deploy_spot_comp, j, i, k, users_comp{j,i}(k,5), angle_3dB{j,i}(k,1), Am ); %azimuth radiation pattern
        A_H=10*log10(A_H); %from linear to dB
        A_V = evalAV(deploy_spot_comp, j, i, k, users_comp{j,i}(k,4), angle_3dB{j,i}(k,2), SLAv ); %elevation radiation pattern
        A_V=10*log10(A_V); %from linear to dB
        A_tx=A_H+A_V; %overall radiation pattern in dB
        A_tx(A_tx<-Am) = -Am; %overall radiation pattern
        %compute Beamforming Gain
        BF = 10*log10(real(sinc(((users_comp{j,i}(k,5)-deploy_spot_comp{j,i}(k,5))/(1.13*angle_3dB{j,i}(k,1))).^2))) +...
            10*log10(real(sinc(((users_comp{j,i}(k,4)-deploy_spot_comp{j,i}(k,4))/(1.13*angle_3dB{j,i}(k,2))).^2))) + Gbf;
        BF=real(BF);
        PL=users_comp{j,i}(k,6); %path loss in dB
        P_u_co = txPowerDBm - PL + A_tx + Gtx + BF; %compute desired signal power 
        P_u_co = 10^(P_u_co/10); %from dB to linear
        
        %compute intra sector interference
        for m=1:length(deploy_spot_comp{j,i})
            if m~=k 
                A_H = evalAH(deploy_spot_comp, j, i, m, users_comp{j,i}(k,5), angle_3dB{j,i}(m,1), Am ); %azimuth radiation pattern
                A_H=10*log10(A_H); %from linear to dB
                A_V = evalAV(deploy_spot_comp, j, i, m, users_comp{j,i}(k,4), angle_3dB{j,i}(m,2), SLAv ); %elevation radiation pattern
                A_V=10*log10(A_V); %from linear to dB
                A_tx=A_H+A_V; %overall radiation pattern in dB
                A_tx(A_tx<-Am) = -Am; %overall radiation pattern
                %compute Beamforming Gain
                BF = 10*log10(sinc(((users_comp{j,i}(k,5)-deploy_spot_comp{j,i}(m,5))/(1.13*angle_3dB{j,i}(m,1))).^2)) +...
                    10*log10(sinc(((users_comp{j,i}(k,4)-deploy_spot_comp{j,i}(m,4))/(1.13*angle_3dB{j,i}(m,2))).^2)) + Gbf;
                BF=real(BF);
                %compute interference
                I = txPowerDBm - PL + A_tx + Gtx + BF; %intra sector interference
                I = 10^(I/10); %from dB to linear 
                I_u_co_so = [I_u_co_so, I];
            end
        end
        
        %compute inter sector interference (same cell)
        for l=1:size(deploy_spot_comp,2)
            if l~=i 
                for m=1:length(deploy_spot_comp{j,l})
                    center=[BS(j,1),BS(j,2)]; %position of current BS
                    steer = evalsteer(l, pos, center); %computes steer angle
                    A_H = evalAH(deploy_spot_comp, j, l, m, steer, angle_3dB{j,l}(m,1), Am ); %azimuth radiation pattern
                    A_H=10*log10(A_H); %from linear to dB
                    A_V = evalAV(deploy_spot_comp, j, l, m, users_comp{j,i}(k,4), angle_3dB{j,l}(m,2), SLAv ); %elevation radiation pattern
                    A_V=10*log10(A_V); %from linear to dB
                    A_tx=A_H+A_V; %overall radiation pattern in dB
                    A_tx(A_tx<-Am) = -Am; %overall radiation pattern
                    %compute Beamforming Gain
                    BF = 10*log10(real(sinc(((steer-deploy_spot_comp{j,l}(m,5))/(1.13*angle_3dB{j,l}(m,1))).^2))) +...
                        10*log10(real(sinc(((users_comp{j,i}(k,4)-deploy_spot_comp{j,l}(m,4))/(1.13*angle_3dB{j,l}(m,2))).^2))) + Gbf;
                    BF=real(BF);
                    %compute interference
                    I = txPowerDBm - PL + A_tx + Gtx + BF; %inter sector interference
                    I = 10^(I/10); %from dB to linear
                    I_u_co_s = [I_u_co_s, I];
                end
            end
        end
        
        %compute inter cell interference
        for n=2:size(deploy_spot_comp,1)
            for l=1:size(deploy_spot_comp,2)
                for m=1:length(deploy_spot_comp{n,l})
                    center=[BS(n,1),BS(n,2)]; %position of current BS
                    d=norm(pos-center); %2D distance between user and BS
                    d_eff=sqrt((effective_h)^2 + (d)^2); %effective distance
                    tilt=atan2d(effective_h,d); %computes tilt angle
                    steer = evalsteer(l, pos, center); %computes steer angle
                    A_H = evalAH(deploy_spot_comp, n, l, m, steer, angle_3dB{n,l}(m,1), Am ); %azimuth radiation pattern
                    A_H=10*log10(A_H); %from linear to dB
                    A_V = evalAV(deploy_spot_comp, n, l, m, tilt, angle_3dB{n,l}(m,2), SLAv ); %elevation radiation pattern
                    A_V=10*log10(A_V); %from linear to dB
                    A_tx=A_H+A_V; %overall radiation pattern in dB
                    A_tx(A_tx<-Am) = -Am; %overall radiation pattern
                    %compute Beamforming Gain
                    BF = 10*log10(real(sinc(((steer-deploy_spot_comp{n,l}(m,5))/(1.13*angle_3dB{n,l}(m,1))).^2))) +...
                        10*log10(real(sinc(((tilt-deploy_spot_comp{n,l}(m,4))/(1.13*angle_3dB{n,l}(m,2))).^2))) + Gbf;
                    BF=real(BF);
                    %compute path loss
                    PL_ext=21*log10(d_eff) + 32.4 + 20*log10(fq);
                    if NLOS_on==1
                        PL_nlos_ext=35.3*log10(d_eff)+22.4+21.3*log10(fq)-0.3*(UE_h-1.5);
                        PL=max(PL_nlos_ext,PL_ext);
                    end
                    %compute interference
                    I = txPowerDBm - PL_ext + A_tx + Gtx + BF; %inter sector interference
                    I = 10^(I/10); %from dB to linear
                    I_u_c_s = [I_u_c_s, I];
                end
            end
        end
        
        %compute SINR 
        if int_intrasec==1
            SINR = P_u_co / (sum(I_u_co_so)+sum(I_u_co_s)+sum(I_u_c_s)+rxNoisePower);
        else
            SINR = P_u_co / (sum(I_u_co_s)+sum(I_u_c_s)+rxNoisePower);
        end
        instant_ach_rate = real(log2(1+(SINR)));
        tp_user{j,i}(k)=bw*instant_ach_rate; %throughput compued over user
        
        % deviation of UE values wrt deploy spot
        delta_tp{j,i}(k)=tp_ds{j,i}(k)-tp_user{j,i}(k); %delta throughput between deploy spot and user
        center=[BS(j,1),BS(j,2)]; %position of current BS
        pos_ds=[deploy_spot_comp{j,i}(k,1),deploy_spot_comp{j,i}(k,2)]; %position of deployment spot
        if norm(center-pos) < norm(center-pos_ds)
            delta_dist{j,i}(k)=-(norm(pos_ds-pos));
        else
            delta_dist{j,i}(k)=norm(pos_ds-pos);
        end
        dist_user{j,i}(k)=norm(center-pos); %2D distance between user and BS
        
    end
end

%% plot

if plot_on==1
    for i=1:length(BS)
        imm1=plot(sector{i,1});
        hold on
        imm2=plot(sector{i,2});
        hold on
        imm3=plot(sector{i,3});
        hold on
        axis equal;
        if fig_style == 1 && i~=1
            imm1.FaceColor = 'none';
            imm2.FaceColor = 'none';
            imm3.FaceColor = 'none';
        end
    end
    for i=1:size(deploy_spot_comp,2)
        for k=1:length(deploy_spot_comp{1,i})
            if fig_style~=1 
            plot(deploy_spot_comp{1,i}(k,1),deploy_spot_comp{1,i}(k,2),'b*');
            hold on;
            end
            plot(users_comp{1,i}(k,1),users_comp{1,i}(k,2),'r.');
            hold on;
        end
    end
end
% toc

%% Angle implementation for Point Source Model

%Compute number of measure spot per each sector
numb_ms=0;
for i=1:size(measure_spot,2)
    numb_ms=numb_ms+length(measure_spot{1,i}); 
end
Seq_ms{1,numb_ms}=[];
c=1;
emf_ms=[];
data_ms=[]; 

%%No Beamforming adopt point-source mode of ITU: gNB is collapsed into an 
%%omnidirectional base station always radiating at txPower
if fixed_emf == 1
    EIRP = txPower*G;
    EIRP_db = 10*log10(EIRP);
    for i=1:size(measure_spot,2)
        for k=1:length(measure_spot{1,i})
            pos=[measure_spot{1,i}(k,1),measure_spot{1,i}(k,2)];
            for j=1:length(BS)
                center=[BS(j,1),BS(j,2)]; %postion of current BS
                dist_ms=norm(pos-center); %2D distance between measure spot and BS
                dist_eff=sqrt((dist_ms)^2+(effective_h)^2); %effective distance
                Seq(j)=(EIRP/(4*pi*(dist_eff)^2)); 
            end
            emf_ms =[emf_ms; [i,measure_spot{1,i}(k,1),measure_spot{1,i}(k,2),sqrt(sum(Seq)*Z0)]];
        end
    end
end
save('fixed_emf_data.mat','emf_ms');

%%EMF computation when Beamforming active
if measure_on==1
    for i=1:size(measure_spot,2)
        for k=1:length(measure_spot{1,i})
            pos=[measure_spot{1,i}(k,1),measure_spot{1,i}(k,2)]; %position of current measure spot
            for j=1:length(BS)
                center=[BS(j,1),BS(j,2)]; %coordinates of current BS
                dist_ms=norm(pos-center); %2D distance between measure spot and BS
                dist_eff=sqrt((dist_ms)^2+(effective_h)^2);%effective distance
                tilt_ms=atan2d(effective_h,dist_2D); %computes tilt angle
                if j==1
                    for sec=1:size(measure_spot,2)
                        steer = evalsteer(sec, pos, center); %compute steer angle
                        for ds=1:length(deploy_spot_comp{j,sec})
                            A_H = evalAH(deploy_spot_comp, j, sec, ds, steer, angle_3dB{j,sec}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, sec, ds, tilt_ms, angle_3dB{j,sec}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end
                end
                if j==2
                    %sector 2
                    steer = evalsteer(2, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,2})
                            A_H = evalAH(deploy_spot_comp, j, 2, ds, steer, angle_3dB{j,2}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 2, ds, tilt_ms, angle_3dB{j,2}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                    end
                elseif j==3
                    %sector 2
                    steer=evalsteer(2, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,2})
                        if deploy_spot_comp{j,2}(ds,3)==0
                           A_H = evalAH(deploy_spot_comp, j, 2, ds, steer, angle_3dB{j,2}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 2, ds, tilt_ms, angle_3dB{j,2}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end
                    %sector 3
                    steer=evalsteer(3, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,3})
                        if deploy_spot_comp{j,2}(ds,3)==1
                            A_H = evalAH(deploy_spot_comp, j, 3, ds, steer, angle_3dB{j,3}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 3, ds, tilt_ms, angle_3dB{j,3}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end
                elseif j==4
                    %sector 3
                    steer=evalsteer(3, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,3})
                           A_H = evalAH(deploy_spot_comp, j, 3, ds, steer, angle_3dB{j,3}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 3, ds, tilt_ms, angle_3dB{j,3}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                    end
                elseif j==5
                    %sector 1
                    steer=evalsteer(1, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,1})
                        if deploy_spot_comp{j,2}(ds,3)==1
                            A_H = evalAH(deploy_spot_comp, j, 1, ds, steer, angle_3dB{j,1}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 1, ds, tilt_ms, angle_3dB{j,1}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end
                    %sector 3
                    steer=evalsteer(3, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,3})
                        if deploy_spot_comp{j,2}(ds,3)==0
                            A_H = evalAH(deploy_spot_comp, j, 3, ds, steer, angle_3dB{j,3}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 3, ds, tilt_ms, angle_3dB{j,3}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end
                elseif j==6
                    %sector 1
                    steer=evalsteer(1, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,1})
                            A_H = evalAH(deploy_spot_comp, j, 1, ds, steer, angle_3dB{j,1}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 1, ds, tilt_ms, angle_3dB{j,1}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                    end
                elseif j==7
                    %sector 1
                    steer=evalsteer(1, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,1})
                        if deploy_spot_comp{j,2}(ds,3)==0
                            A_H = evalAH(deploy_spot_comp, j, 1, ds, steer, angle_3dB{j,1}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 1, ds, tilt_ms, angle_3dB{j,1}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end
                    %sector 2
                    steer=evalsteer(2, pos, center);%compute steer angle
                    for ds=1:length(deploy_spot_comp{j,2})
                        if deploy_spot_comp{j,2}(ds,3)==1
                            A_H = evalAH(deploy_spot_comp, j, 2, ds, steer, angle_3dB{j,2}(1,1), Am );
                            A_V = evalAV(deploy_spot_comp, j, 2, ds, tilt_ms, angle_3dB{j,2}(1,2), SLAv );
                            Seq=((PG_4pi)/(dist_eff)^2)*(A_H)^2*(A_V)^2;
                            Seq_ms{c}=[Seq_ms{c}; Seq];
                        end
                    end

                end
            end
            c=c+1;
            data_ms =[data_ms; [i,measure_spot{1,i}(k,1),measure_spot{1,i}(k,2)]];
        end
    end
    
    Seq_sum=[];
    for i=1:size(Seq_ms,2)
        Seq_sum=[Seq_sum sum(Seq_ms{1,i})];
        data_ms(i,4) =sqrt(sum(Seq_ms{1,i})*Z0); %compute EMF
    end
    
    E=sqrt(Seq_sum*Z0);                         % compute electric field [V/m]
    avg_emf=mean(E,'omitnan');                  % avg
    SEM = std(E,'omitnan')/sqrt(length(E));     % Standard Error
    ts = tinv([0.025  0.975],length(E)-1);      % T-Score
    emf_CI = avg_emf + ts*SEM;                  % Confidence Intervals

    fname=join(['Data\emf',num2str(accuracy),'.mat']);
    save(fname,'avg_emf','emf_CI','E');
    save('angle_3dB', 'angle_3dB');
    save('emf_data.mat','data_ms');
end

toc

%% FUNCTIONS

%%A_H and A_V are computed according to
%%https://ieeexplore.ieee.org/document/7465794
function A_H = evalAH(deploy_spot_comp, cell_id, sector_id, ds, steer, az3dB, Am )
    A_H = 12*(((steer-deploy_spot_comp{cell_id,sector_id}(ds,5))/az3dB).^2);  %Horizontal pattern -da steer
    A_H=-(min(A_H,Am));
    A_H=10^(A_H/10);
end
function A_V = evalAV(deploy_spot_comp, cell_id, sector_id, ds, tilt_ms, el3dB, SLAv )
    A_V = 12*(((tilt_ms-deploy_spot_comp{cell_id,sector_id}(ds,4))/el3dB).^2); %vertical pattern -da tilt
    A_V=-(min(A_V,SLAv)); %saturated at maximum attenuation
    A_V=10^(A_V/10);
end

%evalsteer computes steer angle which is the angle computed w.r.t. sector bisector
%sector 1 -> bisector=90°, sector 2 -> bisector = 210°, sector 3 -> bisector=330°
function steer = evalsteer(sector_id, pos, center)
sign=0;
%Compute the angle of the spot characterized by coordinates in "pos" vector
%in a cartesian coordinate system centered in the base station 
%characterized by coordinates in "center" vector
angle = atan2d(pos(1,2)-center(1,2), pos(1,1)-center(1,1));

%Report the angle in a positive 0-360° notation
if angle<0
    angle= 360+angle;
    sign=1;
end

if sector_id ==1
    if  angle> 270 && sign==1
        angle = angle - 360;
    end
    steer= 90-angle;
elseif sector_id ==2
        if  angle < 30 
        angle = 360 + angle;
        end
    steer = 210-angle; 
elseif sector_id ==3
    if  angle < 150 && sign==0
        angle=360+angle;
    end
    steer = 330-angle;
end
end

