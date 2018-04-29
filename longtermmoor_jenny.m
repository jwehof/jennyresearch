%% Connecticut River Long-term Analysis - Jenny Wehof
% Uses the long-term moored CTD data from the CT River. Option of using a
% low-pass filter to remove semidiurnal tidal variation. Salinity data and
% river discharge from USGS at Thompsonville is used to find the
% relationship between salinity intrusion length and river discharge.
clear all; close all;

% Tidally filter data? Yes, set f=1
f=1;
% Produce time series plots? Yes, set p=1
p=0;

% load the data
load 'C:\Data\CTriv\longterm\mooredCtds.mat'
sitename=char('Deep River','Essex','Old Lyme','Saybrook','Deep River','Essex','Old Lyme','Saybrook','FZ5','FZ5');
% remove bad wl:
wl(wl<-2)=NaN;

% dates as year-day for each year
yrday=NaN(size(dn));
k = dn >= datenum(2012,01,01,00,00,00) & dn < datenum(2013,01,01,00,00,00);
yrday(k) = dn(k)-datenum(2012,01,01,00,00,00);
k = dn >= datenum(2013,01,01,00,00,00) & dn < datenum(2014,01,01,00,00,00);
yrday(k) = dn(k)-datenum(2013,01,01,00,00,00);
k = dn >= datenum(2014,01,01,00,00,00) & dn < datenum(2015,01,01,00,00,00);
yrday(k) = dn(k)-datenum(2014,01,01,00,00,00);
k = dn >= datenum(2015,01,01,00,00,00) & dn < datenum(2016,01,01,00,00,00);
yrday(k) = dn(k)-datenum(2015,01,01,00,00,00);

% convert lat/long coordinates to x-y and cross- and along-channel planes 
[x,y,xc,yc]=ctriv_xy(lat,lon); % also requires 'ct_thalweg.mat'

run defaultfigure.m
yr=[2012,repmat(2013,1,6),repmat(2014,1,6),repmat(2015,1,6)];
mo=[11,1:2:12,1:2:12,1:2:12];
DateLabel=datenum(yr,mo,ones(size(yr))); % datetick labels for the 1st of every other month, Nov 2012-Nov 2015


%% USGS Data
% USGS: ss (salinity at surface), sb (salinity at bottom), wl and ssc at Essex
% Interpolate to be on same time stamp
%datstr='201208to201406';
datstr='200710to201512';
A=load(['C:\Data\CTriv\usgs\usgs_essex_201002to201512.mat']);
wl_e=interp1(A.dnwl,A.wl,dn)';
ss_e=interp1(A.dn,A.ss,dn)';
sb_e=interp1(A.dn,A.sb,dn)';
ssc_e=interp1(A.dn,A.obs,dn)';
tb_e=interp1(A.dn,A.tb,dn)';
ts_e=interp1(A.dn,A.ts,dn)';
lat_e=A.lat;
lon_e=A.lon;

A=load(['C:\Data\CTriv\usgs\usgs_oldlyme_200710to201512.mat']);
wl_ol=interp1(A.dnwl,A.wl,dn)';
ss_ol=interp1(A.dn,A.ss,dn)';
sb_ol=interp1(A.dn,A.sb,dn)';
ssc_ol=interp1(A.dn,A.cb,dn)';
tb_ol=interp1(A.dn,A.tb,dn)';
ts_ol=interp1(A.dn,A.ts,dn)';
lat_ol=A.lat;
lon_ol=A.lon;

% usgs discharge @ thompsonville
thompsonville=load(['C:\Data\CTriv\usgs\usgs_thompsonville_' datstr '.mat']);
qr=interp1(thompsonville.dnwl,thompsonville.qr,dn); % put qr on same timestamp as CTD data

%% Calculate tidal amplitude
% fs=1/300; % sampling frequency, Hz
% wlmin=tmin(wl,fs);
% wlmax=tmax(wl,fs);
% amp=wlmax-wlmin;
% ampfilt=pl33tn(amp,1/12,33);

% save outputs to save time
load amp.mat
load ampfilt.mat
plot(dn,ampfilt(:,2));
title('Tidal Amplitude at Essex')
ylabel('tidal amplitude [m]')
xlabel('Date')
set(gca,'XTick',DateLabel);
datetick('x',3,'keepticks'); grid on;

%% Filter Data
% pl33 Filter from p. 21, Rosenfeld (1983), WHOI Technical Report 85-35
% Uses default filter half-amplitude period of 33 hr to remove
% diurnal/semidiurnal tides
% keep raw data
sraw=s; nturaw=ntu; wlraw=wl; traw=t;
if f==1
    s=pl33tn(s,1/12,33);
    ntu=pl33tn(ntu,1/12,33);
    wl=pl33tn(wl,1/12,33);
    t=pl33tn(t,1/12,33);
end

%% Plots: time series data
% Plots of salinity and turbidity (as ntu) with river discharge and water
% level at Essex

% SALINITY & DISCHARGE
if p==1
    figure
    subplot(3,1,1);
    plot(dn,qr); ylabel('Discharge [m^3 s^-^1]'); ylim([0 2800]);
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,2); 
    plot(dn,s(:,isbot==0)), ylabel('Surface Salinity [psu]'); ylim([0 35]);
    legend(staname(isbot==0),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,3); 
    plot(dn,s(:,isbot==1)); xlabel('Month'); ylabel('Bottom Salinity [psu]'); ylim([0 35]);
    legend(staname(isbot==1),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    % TURBIDITY & DISCHARGE
    figure
    subplot(3,1,1);
    plot(dn,qr);ylabel('Discharge, [m^3 s^-^1]');ylim([0 2800]);
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,2); 
    plot(dn,ntu(:,isbot==0)); ylabel('Surface Turbidity [ntu]'); ylim([0 450]);
    legend(staname(isbot==0),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,3); 
    plot(dn,ntu(:,isbot==1)); xlabel('Month'); ylabel('Bottom Turbidity [ntu]'); ylim([0 450]);
    legend(staname(isbot==1),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    %SALINITY & WATER LEVEL
    figure
    subplot(3,1,1);
    plot(dn,ampfilt(:,8));ylabel('Tidal Range [m]');
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,2); 
    plot(dn,s(:,isbot==0)); ylabel('Surface Salinity [psu]');ylim([0 35]);
    legend(staname(isbot==0),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,3); 
    plot(dn,s(:,isbot==1));xlabel('Month'); ylabel('Bottom Salinity [psu]');ylim([0 35]);
    legend(staname(isbot==1),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    %TURBIDITY & WATER LEVEL
    figure
    subplot(3,1,1);
    plot(dn,wl);ylabel('Water Level [m]');
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,2); 
    plot(dn,ntu(:,isbot==0)); ylabel('Surface Turbidity [ntu]');ylim([0 450]);
    legend(staname(isbot==0),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,3); 
    plot(dn,ntu(:,isbot==1)); xlabel('Month'); ylabel('Bottom Turbidity [ntu]'); ylim([0 450]);
    legend(staname(isbot==1),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;
end

%% compare salinity at OL to USGS data - look for fouling in saybrook

figure
% subplot(2,1,1)
plot(dn,sb_ol,'k');
hold on;
plot(dn,sraw(:,7),'b')
plot(dn,sraw(:,8),'r')
% plot(dn,ss_ol,'--k');
% plot(dn,sraw(:,3),'--b');hold on
legend('USGS bot','Old Lyme bot','Saybrook bot')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('s [psu]')
title('Salinity at Old Lyme')

% difference
subplot(2,1,2)
plot(dn,sb_ol-sraw(:,7),'r'); hold on
% plot(dn,ss_ol-sraw(:,3),'b'); 
title('USGS salinity - WHOI salinity, Old Lyme')
legend('bottom','surface')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('(s_U_S_G_S - s_W_H_O_I')

% filtered
figure
plot(dn,pl33tn(sb_ol,1/12,33),'k');hold on
plot(dn,s(:,7),'b');
plot(dn,s(:,8),'r');
% plot(dn,s(:,3),'--b');hold on
% plot(dn,pl33tn(ss_ol,1/12,33),'--k');
legend('USGS bot','Old Lyme bot','Saybrook bot','USGS surf','Old Lyme surf')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('s [psu]')
title('Tidally Filtered Bottom Salinity')

% filtered difference
figure
plot(dn,pl33tn(sb_ol,1/12,33)-s(:,7),'r'); hold on
plot(dn,pl33tn(ss_ol,1/12,33)-s(:,3),'b'); 
title('USGS salinity - Our salinity, Old Lyme')
legend('bottom','surface')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('(s_U_S_G_S - s_W_H_O_I)')
title('tidally filtered salinity difference between USGS and WHOI')

% filtered surface sensor at OL
figure
% subplot(2,1,1)
plot(dn,pl33tn(ss_ol,1/12,33),'k');hold on;
plot(dn,s(:,3),'b');
plot(dn,s(:,4),'r');
legend('USGS surf','Old Lyme surf','Saybrook surf')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('s [psu]')
title('Tidally Filtered Salinity')

% difference
% subplot(2,1,2)
plot(dn,pl33tn(ss_ol,1/12,33)-s(:,3),'r'); hold on
plot(dn,pl33tn(ss_ol,1/12,33)-s(:,4),'b'); 
title('USGS salinity - WHOI salinity')
legend('Old Lyme','Saybrook')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('(s_U_S_G_S - s_W_H_O_I')


% water level
figure
plot(dn,pl33tn(wl_ol,1/12,33)-wl(:,7),'r'); hold on
title('USGS wl - Our wl, Old Lyme')
datetick('x',12)

% temperature time series
figure
plot(dn,traw(:,8)); hold on
plot(dn,tb_ol);
title('Bottom Temperature at Old Lyme')
legend('WHOI','USGS')
% title('Temperature at Saybrook')
ylabel('degrees Celcius')
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')

%% T-S relationship - tidal time scale
% Old Lyme
figure
subplot(2,2,1)
date = datenum(2015,2,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Old Lyme,' datestr(date)])
subplot(2,2,2)
date = datenum(2015,5,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Old Lyme,' datestr(date)])
subplot(2,2,3)
date = datenum(2015,8,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Old Lyme,' datestr(date)])
subplot(2,2,4)
date = datenum(2015,11,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Old Lyme,' datestr(date)])

% Saybrook
figure
subplot(2,2,1)
date = datenum(2014,1,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(traw(dn>=tidetime(1) & dn<=tidetime(2),8),sraw(dn>=tidetime(1) & dn<=tidetime(2),8))
hold on;
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Saybrook,' datestr(date)])
subplot(2,2,2)
date = datenum(2014,4,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(traw(dn>=tidetime(1) & dn<=tidetime(2),8),sraw(dn>=tidetime(1) & dn<=tidetime(2),8))
hold on;
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Saybrook,' datestr(date)])
subplot(2,2,3)
date = datenum(2014,7,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(traw(dn>=tidetime(1) & dn<=tidetime(2),8),sraw(dn>=tidetime(1) & dn<=tidetime(2),8))
hold on;
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Saybrook,' datestr(date)])
subplot(2,2,4)
date = datenum(2014,10,1);
tidetime = [date-datenum(0,0,0,6.21,0,0) date+datenum(0,0,0,6.21,0,0)];
scatter(traw(dn>=tidetime(1) & dn<=tidetime(2),8),sraw(dn>=tidetime(1) & dn<=tidetime(2),8))
hold on;
scatter(tb_ol(dn>=tidetime(1) & dn<=tidetime(2)),sb_ol(dn>=tidetime(1) & dn<=tidetime(2)))
ylabel('Salinity [psu]'); xlabel('Temperature [degC]');
title(['T-S Plot at Saybrook,' datestr(date)])
legend('Saybrook','Old Lyme (USGS)')

% find T-S for all data, regression and r-sq to see if we can reconstruct s
j=find(dn<=dn(1)+datenum(0,0,0,6.21,0,0));
TSreg = NaN(length(dn),2) ; TSrsq = NaN(size(dn));
for i=j(end)+1:length(dn)-j-1
    tidetime = [dn(i)-datenum(0,0,0,6.21,0,0) dn(i)+datenum(0,0,0,6.21,0,0)];
    k = dn>=tidetime(1) & dn<=tidetime(2);
    T = tb_ol(k); S = sb_ol(k);
    [p1,s1,mu1] = polyfit(T,S,1);
    TSreg(i,:) = p;
    TScorr = corrcoef(T,S); TSrsq(i) = TScorr(2)^2;
end
    


%% Data Preparation
% Prepare the bottom salinity data to be used in the calculation of the
%   salinity intrusion length. Fills in some holes in the data or corrects negative
%   values produced from filtering. Plots salinity contours with discharge to see trend. 'sbot' stations are in columns, 1. Deep
%   River, 2. Essex, 3. Old Lyme, 4. Saybrook, 5. Frontal Zone 5

sbot=s(:,isbot==1);
n=isnan(sbot);
sbot(n(:,1)==1,1)=s(n(:,1)==1,1); % for missing Deep River bottom data, use Deep River surface (well-mixed)
sbot(sbot < 0)=0; % set all negative salinities produced in filtering to 0

% salinities at FZ5 which are less than Deep River get set to the Deep
% River value
k=find(sbot(:,5) <= sbot(:,1)); 
sbot(k,5)=sbot(k,1);
% salinities at Essex which are less than FZ5 get set to FZ5 value
k=find(sbot(:,2) <= sbot(:,5)); 
sbot(k,2)=sbot(k,5);
% salinities at Essex which are less than Deep River get set to DR value
k=find(sbot(:,2) <= sbot(:,1)); 
sbot(k,2)=sbot(k,1);

% distance of each bottom ctd from the mouth, in km
xbot=xc(isbot==1)/1000;
xbot=repmat(xbot,length(dn),1);

%% stratification

ds=sbot;
ds(:,1)=sbot(:,1)-s(:,1); % Deep River
ds(:,2)=sbot(:,2)-s(:,2); % Essex
ds(:,3)=sbot(:,3)-s(:,3); % Old Lyme
ds(:,4)=sbot(:,4)-s(:,4); % Saybrook
ds(:,5)=sbot(:,5)-s(:,10); % FZ5

dstmax=tmax(ds,fs);
dstmin=tmin(ds,fs);
dstrange=dstmax-dstmin;

% plot stratification with discharge and tidal range
% all data together
figure
subplot(3,1,1);
plot(dn,qr,'k','LineWidth',1.5); ylabel('Discharge [m^3 s^-^1]','FontSize',12); ylim([0 2800]);
set(gca,'XTick',DateLabel);
datetick('x',3,'keepticks'); grid on;
title('Tidal Minimum Stratification','FontSize',16);

subplot(3,1,2);
plot(dn,ampfilt(:,8), 'k','LineWidth',1.5); ylabel('Tidal Range [m]','FontSize',12);hold on;
set(gca,'XTick',DateLabel);
datetick('x',3,'keepticks'); grid on;

subplot(3,1,3); 
plot(dn,pl33tn(dstmin(:,2:end),1/12,33)); ylabel('\Delta s = s_b_o_t - s_s_u_r_f','FontSize',12);ylim([0 35]);
legend('Essex','Old Lyme','Saybrook','FZ5','orientation','horiz'); legend boxoff;
set(gca,'XTick',DateLabel);
datetick('x',3,'keepticks'); grid on;

% 1 year at a time
year=2013; k = dn >= datenum(2013,01,01,00,00,00) & dn < datenum(2014,01,01,00,00,00); DateLabel2=[datenum(repmat(2013,1,12),1:12,ones(1,12)) datenum(2014,1,1)];
% year=2014; k = dn >= datenum(2014,01,01,00,00,00) & dn < datenum(2015,01,01,00,00,00); DateLabel2=[datenum(repmat(2014,1,12),1:12,ones(1,12)) datenum(2015,1,1)];
% year=2015; k = dn >= datenum(2015,01,01,00,00,00) & dn < datenum(2016,01,01,00,00,00); DateLabel2=[datenum(repmat(2015,1,7),1:7,ones(1,7))];

figure
subplot(3,1,1);
plot(dn(k),qr(k),'k','LineWidth',1.5); ylabel('Discharge [m^3 s^-^1]','FontSize',12); ylim([0 2800]);
set(gca,'XTick',DateLabel2);
datetick('x',3,'keepticks'); grid on;
title('Tidal Minimum Stratification','FontSize',16);

subplot(3,1,2);
plot(dn(k),ampfilt(k,8), 'k','LineWidth',1.5); ylabel('Tidal Range [m]','FontSize',12);hold on;
set(gca,'XTick',DateLabel2);
datetick('x',3,'keepticks'); grid on;

subplot(3,1,3); 
plot(dn(k),pl33tn(dstmin(k,2:end),1/12,33)); ylabel('\Delta s = s_b_o_t - s_s_u_r_f','FontSize',12);ylim([0 35]);
legend('Essex','Old Lyme','Saybrook','FZ5','orientation','horiz'); legend boxoff;
set(gca,'XTick',DateLabel2);
datetick('x',3,'keepticks'); grid on;
xlabel(year,'FontSize',14,'FontWeight','bold');



%% plot salinity contour time series
p=0;
if p==1
    dncont=[dn' dn' dn' dn' dn'];
    scont=[0 5 10 15 20 25 30];
    figure
    subplot(2,1,1)
    [C,h] = contourf(dncont,[20*ones(length(xbot),1) xbot(:,1:4)],[zeros(length(sbot),1) sbot(:,1:4)],scont,'LineStyle','none');
    ylabel('Distance from mouth [km]','FontSize',14); 
    hold on;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;
    plot(repmat(dn(1),1,4),xbot(1,1:4),'r+')
    xlim([datenum(2013,05,01,0,0,0) datenum(2015,06,01,0,0,0)])
    h=colorbar; ylabel(h,'psu');

    subplot(2,1,2)
    plot(dn,qr,'k'); 
    ylabel('Discharge [m^3 s^-^1]','FontSize',14);
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;
    xlim([datenum(2013,05,01,0,0,0) datenum(2015,06,01,0,0,0)])
    colorbar
    
% to use this, change total subplots previously to 3    
%     subplot(3,1,3)
%     [ax,h1,h2] = plotyy(dn,ntu(:,6),dn,s(:,6));
%     set(gca,'XTick',DateLabel);
%     datetick('x',12,'keepticks'); grid on;
%     xlim([datenum(2013,05,01,0,0,0) datenum(2014,06,01,0,0,0)])
%     xlabel('Date','FontSize',14); 
%     set(ax(1),'YLim',[0 200],'Box','off','YTick',[0:50:200])
%     set(ax(2),'YLim',[0 40],'YTick',[0:10:40])
%     colorbar
end

%% Salinity Intrusion
% Calculate the length of salinity intrusion using salinity data from moored
%   stations.  Only use interpolation if the salinity at Old Lyme is 2 psu or greater;
%   otherwise it will not work and returns X2 values less than zero since exponential only can work for the tail end of the data. 

% Exponential for tail end of data - to get X2
lns=log(sbot);
n=~isnan(lns);

c=NaN(size([dn' dn'])); % polyfit coefficients
rsq_s_X2=NaN(size(dn')); % r-squared for each polyfit
for i=1:length(dn');
    if sbot(i,3) >= 2.0; 
        c(i,:)=polyfit(lns(i,n(i,:)),xbot(i,n(i,:)),1); % only uses non-NaN values
        rsq=corrcoef(lns(i,n(i,:)),xbot(i,n(i,:)));
        rsq_s_X2(i)=rsq(1,2)^2;
    end
end

% use polyfit coefficients to get X2
S=2;
for i=1:length(dn')
    X2(i,1)=polyval(c(i,:),log(S));
end

% bad data - when X2 is outside the range of xbot values
X2(X2>20 | X2<0)=NaN;
figure

% plot X2 with the salinities
subplot(2,1,1);
plot(dn,X2,'LineWidth',1.1); ylabel('Salinity Intrusion, X_2');
set(gca,'XTick',DateLabel);
datetick('x',3,'keepticks'); grid on;
subplot(2,1,2);
plot(dn,s(:,isbot==1)); ylabel('Bottom Salinity [psu]'); ylim([0 35]);
legend(staname(isbot==1),'orientation','horiz'); legend boxoff;
set(gca,'XTick',DateLabel);
datetick('x',3,'keepticks'); grid on;% 

% plot discharge with X2 and the r-sq
figure
% subplot(3,1,1);
% plot(dn,qr,'linewidth',1.5);ylabel('Discharge [m^3 s^-^1]');datetick('x',3);grid on;
% subplot(3,1,2); plot(dn,X2); grid on; datetick('x',3);
% ylabel('X_2 [km]'); 
% subplot(3,1,3); 
scatter(dn,rsq_s_X2,'.'); grid on; set(gca,'XTick',DateLabel); datetick('x',3,'keepticks');
ylabel('R^2'); 


%% Relationship between X2 and Qr

% % Linear fit
% for idx=1:289
%     qrlag=qr(1:(end-(idx-1)))';X2lag=X2(idx:end);
%     m=~isnan(X2lag); % exclude when X2 is NaN
%     qrlag=qrlag(m); X2lag=X2lag(m);
%     k=qrlag <= 1000;
%     X2lag=X2lag(k); qrlag=qrlag(k); % exclude when discharge is greater than 1000 cms
%     reg(idx,:)=polyfit(qrlag,X2lag,1);
%     rsq=corrcoef(qrlag,X2lag);
%     rsq_qr_X2(idx,1)=rsq(1,2)^2;
% end
% disp('rsq - linear =')
% disp(max(rsq_qr_X2))
% 
% % Plot the R-sq vs lag time to find best fit
% lagtime=linspace(0,4,289)';
% figure
% plot(lagtime,rsq_qr_X2);xlabel('Lag time, days');ylabel('R-square');title('Lagged Discharge Effects on Correlation of X2 and Qr');
% 
% k=find(rsq_qr_X2==max(rsq_qr_X2));
% figure
% hold on; grid on;title('Linear Reg. with lag, Q<1000 cms')
% scatter(qrlag,X2lag,'.','MarkerEdgeColor',[0.5 0.5 0.5]);  plot(qrlag,polyval(reg(k,:),qrlag),'k');
% xlabel('Flow, cms'); ylabel('X2, km')

% Using log/log scale and linear reg
logX2=log10(X2);logqr=log10(qr);
% lagged regression
for idx=1:289
    qrlag=logqr(1:(end-(idx-1)))';X2lag=logX2(idx:end);
    m=~isnan(X2lag); % exclude when X2 is NaN
    qrlag=qrlag(m); X2lag=X2lag(m);
    m=~isnan(qrlag);
    qrlag=qrlag(m); X2lag=X2lag(m);
    k=qrlag <= log10(1000);
    %X2lag=X2lag(k); qrlag=qrlag(k); % exclude when discharge is greater than 1000 cms
    reglog(idx,:)=polyfit(qrlag,X2lag,1);
    rsqlog=corrcoef(qrlag,X2lag);
    rsq_qr_X2(idx,1)=rsqlog(1,2)^2;
end
disp('rsq-log =')
disp(max(rsq_qr_X2))

% Plot the R-sq vs lag time to find best fit
lagtime=linspace(0,4,289)';
figure
plot(lagtime,rsq_qr_X2);xlabel('Lag time, days');ylabel('R-square');title('Lagged Discharge Effects on Correlation of X2 and Qr');

k=find(rsq_qr_X2==max(rsq_qr_X2));
disp('regression coefficients - log')
reglog(k,:)
figure
hold on; grid on;title('log-log Linear Reg. with lag, tidally filtered data')
scatter(qrlag,X2lag,'.','MarkerEdgeColor',[0.5 0.5 0.5]);  plot(qrlag,polyval(reglog(k,:),qrlag),'k');
xlabel('log of Flow, cms'); ylabel('log of X2, km')






% qrbin=0:50:2700;
% X2binavg(1,1)=NaN;
% for ii=2:length(qrbin)
%     X2binavg(1,ii)=nanmean(X2(qr>qrbin(ii-1) & qr<=qrbin(ii)));
% end
% figure
% loglog(qrbin,X2binavg,'db')
% xlabel('Q_r');ylabel('X_2')
% c=polyfit(log(qrbin(~isnan(X2binavg))),log(X2binavg(~isnan(X2binavg))),1);
% n=c(1)
% loglog(qrbin,polyval(c,log(qrbin)),'k')

%% Tidal 
% tidal range and tidal average/max/min
% *****CAN'T USE FILTERED DATA!!!******

fs=1/600; % sampling frequency
wlmin=tmin(wl,fs);
wlmax=tmax(wl,fs);
amp=wlmax-wlmin;
X2max=tmax(X2,fs);
smax=tmax(s,fs);
ntumin=tmin(ntu,fs);
ntumax=tmax(ntu,fs);
ntumean=tmean(ntu,fs);

% d_eta/dt
dt=(dn(2)-dn(1))*86400; % time between samples, in seconds
d_eta=NaN(size(wl));
for ii=2:length(wl(:,8))
    d_eta(1,:)=NaN;
    d_eta(ii,:)=wl(ii,:)-wl(ii-1,:);
end
ut=d_eta/dt;
plot(dn,ut(:,8));

% PLOTS - X2, Q, TIDAL RANGE
p=1;
if p==1
    ampbin=0.3:0.1:1.6;
    qrbin=0:100:2700;
    for ii=2:length(qrbin)
        for jj=2:length(ampbin)
            X2binavg(1,1)=NaN;
            X2binavg(ii,jj)=nanmean(X2max(amp(qr>qrbin(ii-1) & qr<qrbin(ii),6)>ampbin(jj-1) & amp(qr>qrbin(ii-1) & qr<qrbin(ii),6)<ampbin(jj),1));
        end
    end
    figure
    for ii=2:length(qrbin)
        scatter((repmat(qrbin(ii),length(ampbin)-1,1)),ampbin(2:end),60,log10(X2binavg(ii,2:end)),'fill');
        hold on;
    end
    xlabel('Qr');ylabel('tidal range');
    h=colorbar;
    set(get(h,'ylabel'),'String', 'x2', 'Rotation', 270,'VerticalAlignment', 'Bottom','FontSize',14)
end

k=find(qr>100 & qr<=400);
scatter(amp(k),smax(k,7),'.');
xlabel('tidal range [m]');ylabel('tidal max s [km]');


% Plots - salinity
p=0;
if p==1
    
    figure
    scatter(amp(:,6),smax(:,6),30,qr);xlabel('tidal amplitude, m'); ylabel('tidal max. salinity, psu');
    h= colorbar;
    set(get(h,'ylabel'),'String', 'River Discharge, cms', 'Rotation', 270,'VerticalAlignment', 'Bottom')

    figure
    scatter(X5,smax(:,6),30,qr);xlabel('Distance to 5 psu');ylabel('tidal max. salinity, psu');
    colorbar;

    t=find(qr>400&qr<500);
    figure
    scatter(amp(t,7),X2(t),30,qr(t));ylabel('salinity intrusion (to 2 psu)');xlabel('tidal amplitude, m')
    colorbar;

    %Tidal max salinity
    figure
    subplot(3,1,1);
    plot(dn,qr); ylabel('Discharge, cms'); ylim([0 2500]);
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,2); 
    plot(dn,smax(:,isbot==0)), ylabel('Surface Salinity, psu'); ylim([0 40]);
    legend(staname(isbot==0),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;

    subplot(3,1,3); 
    plot(dn,smax(:,isbot==1)); xlabel('Date'); ylabel('Bottom Salinity, psu'); ylim([0 40]);
    legend(staname(isbot==1),'orientation','horiz'); legend boxoff;
    set(gca,'XTick',DateLabel);
    datetick('x',3,'keepticks'); grid on;
end


% Plots - turbidity

p=0;
if p==1
    k=find(qr>200 & qr<600);
    figure
    scatter(amp(k,7),ntumin(k,7),30,qr(k));
    xlabel('tidal amplitude, m'); 
    ylabel('tidal avg turbidity, ntu');
    title('Old Lyme Bottom Turbidity, 200<Q<600');
    h = colorbar;
    set(get(h,'ylabel'),'String', 'River Discharge, cms', 'Rotation', 270,'VerticalAlignment', 'Bottom')

    figure
    scatter(qr,ntumin(:,2)); 
    xlabel('River Discharge, cms')
    ylabel('Tidal Minimum Surface Turbidity, ntu')
    title('Essex')
    
    
    % bin-average by discharge
    qrbin=0:100:2700;
    clear ntubinavg;
    for ii=2:length(qrbin)
        ntubinavg(1,:)=NaN;
        for jj=1:10
            ntubinavg(ii,jj)=nanmean(ntumin(qr>qrbin(ii-1) & qr<qrbin(ii),jj));
        end
    end
    xtiled=repmat(xc,length(qrbin),1);
    figure
    for ii=5:9; % pick stations to use
        scatter(qrbin,xtiled(:,ii)/1000,100,log10(ntubinavg(:,ii)),'fill')
        text(50,xtiled(1,ii)/1000+0.25,sitename(ii,:),'VerticalAlignment','bottom','FontSize',12);
        hold on
    end
    xlabel('River Discharge at Thompsonville, cms');ylabel('Distance from mouth, km');
    title('Sediment Map using the tidal minimum bottom turbidity')
    cticks=[0.1 1 10 100 300];
    caxis([-1 2]);
    h=colorbar('YTick',log10(cticks),'YTickLabel',cticks);
    set(get(h,'ylabel'),'String', 'OBS Bottom Turbidity, ntu', 'Rotation', 270,'VerticalAlignment', 'Bottom','FontSize',14)
    
    % bin-average by tidal amplitude
    ampbin=0:0.05:2;
    clear ntubinavg;
    for jj=1:10
        ntubinavg(1,:)=NaN;
        for ii=2:length(ampbin)
            ntubinavg(ii,jj)=nanmean(ntumean(amp(:,jj)>ampbin(ii-1) & amp(:,jj)<ampbin(ii),jj));
        end
    end
    xtiled=repmat(xc,length(ampbin),1);
    figure
    for ii=5:9; % pick stations to use
        scatter(ampbin,xtiled(:,ii)/1000,60,log10(ntubinavg(:,ii)),'fill')
        text(0.1,xtiled(1,ii)/1000+0.25,sitename(ii,:),'VerticalAlignment','bottom','FontSize',12);
        hold on
    end
    xlabel('Tidal Amplitude, m');ylabel('Distance from mouth, km');
    title('Sediment Map using the tidal average bottom turbidity')
    cticks=[1 10 100 300];
    caxis([0 2.5]);
    h=colorbar('YTick',log10(cticks),'YTickLabel',cticks);
    colorlabel('OBS Bottom Turbidity, ntu');
end

%% split up flood vs ebb - TURBIDITY
% discharge limits
qlim=[0 250];

% EBB
ampbin=0:0.1:2;
clear ntubinavg;
for jj=1:10
    ntubinavg(1,:)=NaN;
    l=find(qr(:,jj) > qlim(1) & qr(:,jj) < qlim(2));
    for ii=2:length(ampbin)
        k=find(ut(l,jj)<0);
        ntubinavg(ii,jj)=nanmean(ntu(amp(k,jj)>ampbin(ii-1) & amp(k,jj)<ampbin(ii),jj));
    end
end
xtiled=repmat(xc,length(ampbin),1);
figure
for ii=5:9; % pick stations to use
    scatter(ampbin,xtiled(:,ii)/1000,100,log10(ntubinavg(:,ii)),'fill')
    text(0.1,xtiled(1,ii)/1000+0.25,sitename(ii,:),'VerticalAlignment','bottom','FontSize',12);
    hold on
end
xlabel('Tidal Amplitude, m');ylabel('Distance from mouth, km');
title('Bottom turbidity, EBB')
cticks=[1 10 100 300];
caxis([0 2.5]);
h=colorbar('YTick',log10(cticks),'YTickLabel',cticks);
colorlabel('OBS Bottom Turbidity, ntu');

% FLOOD
ampbin=0:0.1:2;
clear ntubinavg;
for jj=1:10
    ntubinavg(1,:)=NaN;
    for ii=2:length(ampbin)
        k=find(ut(l,jj)>0);
        ntubinavg(ii,jj)=nanmean(ntu(amp(k,jj)>ampbin(ii-1) & amp(k,jj)<ampbin(ii),jj));
    end
end
xtiled=repmat(xc,length(ampbin),1);
figure
for ii=5:9; % pick stations to use
    scatter(ampbin,xtiled(:,ii)/1000,100,log10(ntubinavg(:,ii)),'fill')
    text(0.1,xtiled(1,ii)/1000+0.25,sitename(ii,:),'VerticalAlignment','bottom','FontSize',12);
    hold on
end
xlabel('Tidal Amplitude, m');ylabel('Distance from mouth, km');
title('Bottom turbidity, FLOOD')
cticks=[1 10 100 300];
caxis([0 2.5]);
h=colorbar('YTick',log10(cticks),'YTickLabel',cticks);
colorlabel('OBS Bottom Turbidity, ntu');

% ALL
ampbin=0:0.1:2;
clear ntubinavg;
for jj=1:10
    ntubinavg(1,:)=NaN;
    for ii=2:length(ampbin)
        k=find(ut(l,jj)>0 | ut(l,jj)<0);
        ntubinavg(ii,jj)=nanmean(ntu(amp(k,jj)>ampbin(ii-1) & amp(k,jj)<ampbin(ii),jj));
    end
end
xtiled=repmat(xc,length(ampbin),1);
figure
for ii=5:9; % pick stations to use
    scatter(ampbin,xtiled(:,ii)/1000,100,log10(ntubinavg(:,ii)),'fill')
    text(0.1,xtiled(1,ii)/1000+0.25,sitename(ii,:),'VerticalAlignment','bottom','FontSize',12);
    hold on
end
xlabel('Tidal Amplitude, m');ylabel('Distance from mouth, km');
title('Bottom turbidity, ALL')
cticks=[1 10 100 300];
caxis([0 2.5]);
h=colorbar('YTick',log10(cticks),'YTickLabel',cticks);
colorlabel('OBS Bottom Turbidity, ntu');
%% Bottom Turb. vs salinity
scatter(s(k,6),ntu(k,6))
% ylim([0 1500]);
xlabel('salinity'),ylabel('turbidity');
title('Essex Bottom')
%% ADCP data
load('C:\Data\CTriv\moor\adcpSaybrook.mat');

yrdaylag=yrday+(1.5/24);
hold on
plot(yrday,nanmean(ura,2),'k')
xlabel('day')
ylabel('depth-avg velocity')
legend('calculated u=d eta/dt','calculated u with phase shift','ADCP at Saybrook') 





