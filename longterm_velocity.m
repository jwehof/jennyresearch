%% ADCP data
clear
close all

run defaultfigure.m
yr=[2012,repmat(2013,1,6),repmat(2014,1,6),repmat(2015,1,6)];
mo=[11,1:2:12,1:2:12,1:2:12];
DateLabel=datenum(yr,mo,ones(size(yr))); % datetick labels for the 1st of every other month, Nov 2012-Nov 2015


load('C:\Data\CTriv\longterm\adcpSaybrook.mat');

%% retrieve NOAA tide data at New London CT (ID# 8461490) - 
% % 6-min data can only be retrieved in chunks of 31 days, we need Nov 2012 to Dec 2015
% begin_date = string([20121101,20121201,20130101:100:20131201,20140101:100:20141201,20150101:100:20151201]);
% end_date = string([20121130,20121231,20130131,20130228,20130331,20130430,20130531,20130630,20130731,20130831,20130930,20131031,20131130,20131231,20140131,20140228,20140331,20140430,20140531,20140630,20140731,20140831,20140930,20141031,20141130,20141231,20150131,20150228,20150331,20150430,20150531,20150630,20150731,20150831,20150930,20151031,20151130,20151231]);
% wl_noaa = [];
% for i=1:length(begin_date);
%     data=retrieve_noaacoops(begin_date(i),end_date(i),'8461490','predictions','MLLW','metric','lst','csv');
%     data2=retrieve_noaacoops(begin_date(i),end_date(i),'8461490','water_level','MLLW','metric','lst','csv');
%     alldata=[datenum(data.DateTime(1:length(data2.WaterLevel))) data.Prediction(1:length(data2.WaterLevel)) data2.WaterLevel];
%     wl_noaa= [wl_noaa; alldata];
% end
% save wl_noaa

load wl_noaa
%% Calculate tidal amplitude
fs=1/360; % sampling frequency, Hz
wlmin=tmin(wl_noaa(:,2),fs);
wlmax=tmax(wl_noaa(:,2),fs);
amp=wlmax-wlmin;
ampfilt=pl33tn(amp,1/12,60);
dn_noaa=wl_noaa(:,1);
plot(wl_noaa(:,1),ampfilt); hold on
set(gca,'XTick',DateLabel);datetick('x',12,'keepticks')
ylabel('Tidal Amplitude [m]')
title('Spring/Neap cycle at New London')
[spring,sloc]=findpeaks(ampfilt);
[neap,nloc]=findpeaks(-ampfilt);
scatter(wl_noaa(sloc,1),ampfilt(sloc),'k','filled')
scatter(wl_noaa(nloc,1),ampfilt(nloc),'k','filled')

% create sigma layers
nsigma=[1:25]';
sigmah=repmat((da./length(nsigma)),length(nsigma),1);
sigma=sigmah.*repmat(nsigma,1,length(da));
% interpolate velocity onto sigma layers
for ii=1:length(dna)
    ursigma(:,ii)=interp1(za,ura(:,ii),sigma(:,ii));
end


zavgu=nanmean(ura);
zavgufilt=pl33tn(zavgu,1/12,33); % depth-avg, tidally filtered velocity

plot(dna,nanmean(ura),'Color',[0.8 0.8 0.8])
hold on
plot(dna,pl33tn(nanmean(ura),1/12,33),'k','LineWidth',2)
ylabel('depth-avg velocity [m/s]','FontSize',14)
datetick('x',12,'keepticks')

% velocity profiles
% time averaged over deployment
avgu=nanmean(ura,2);
figure
plot([0 0],[0 8],'--k'); hold on
plot(avgu(1:22),za(1:22),'Color',[0.5 0.5 0.5]); 
xlabel('velocity [m/s]','FontSize',14)
ylabel('depth [mab]'); 

% when max vel out
k=find(zavgu==max(zavgu));
plot(ura(:,k),za,'r')

k=find(zavgu==min(zavgu));
plot(ura(:,k),za,'b')
legend('u=0','average profile','max flood','max ebb')