clc;
clear all;


Uhsas_cn_hour_all=[];
Uhsas_sd_hour_all=[];
Uhsas_time_hour=[];
List =dir('E:\wustl\size distribution\target clean\outlier\data_merge\uhsas*.mat');
k =length(List);
for i=1:k
    file_name=List(i).name;
    load(file_name);
    Uhsas_cn_hour_all=[Uhsas_cn_hour_all;uhsas_cn_hour_all(1:end-1,:)];
    Uhsas_sd_hour_all=[Uhsas_sd_hour_all;uhsas_sd_hour_all(1:end-1,:)];
    Uhsas_time_hour=[Uhsas_time_hour;uhsas_time_hour(1:end-1,:)];
end

save('uhsas_clean_data_hourly_2014_2023.mat','Dp_','Uhsas_cn_hour_all', ...
    'Uhsas_sd_hour_all','Uhsas_time_hour');

%% visualization of clean data
clear;

load('uhsas_clean_data_hourly_2014_2023.mat');

daterange = [datenum(2014,1,1,0,0,0) datenum(2024,1,1,0,0,0)];


title_string = {'Size distribution','CN'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1800 700])


% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.3]);

time_uhsas_heated = datenum(Uhsas_time_hour)
idx_uhsas_heated = find(time_uhsas_heated>=daterange(1) & time_uhsas_heated<daterange(2));

x_all = time_uhsas_heated(idx_uhsas_heated);
y1 = Uhsas_sd_hour_all(idx_uhsas_heated,:);

% Data correction
y1(:,52)=(y1(:,51)+y1(:,53))/2;
y1(:,27)=(y1(:,28)+y1(:,26))/2;

PC1 = pcolor(x_all,Dp_,y1');

set(PC1,'EdgeColor','none')
caxis([0 1200])
ax1.YScale = 'log';
set(ax1,'FontSize',12)
xlim(daterange)
datetick('x','yyyy-mm-dd')
% ax1.XTickLabel = '';
ax1.YLabel.String = 'D_p (nm)';
ax1.YLabel.FontSize = 15;
ax1.XAxis.MinorTick = 'off';
ax1.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax1.XAxis.TickDirection = 'out';
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.TickDirection = 'out';
ax1.YLim = [60 1000];
title(title_string{1},'FontSize',15)

h1 = colorbar('v');
h1.Limits = [0 1200];
colormap(jet)
h1.Position = [0.89 0.08 0.015 0.3];
h1.FontSize = 12;
h1.Label.String = 'd{\itN}/d{\itlogD_p} (cm^{–3})';
h1.Label.FontSize = 15;
h1.Ticks = [0:300:1200];
h1.TickDirection = 'out';
h1.TickLength = 0.03;

% 'CN'
ax2 = axes('Position',[0.08 0.5 0.80 0.3]);

y2 = Uhsas_cn_hour_all(idx_uhsas_heated,:);
plot(x_all,y2);
hold on

set(ax2,'FontSize',12)
datetick('x')
xlim(daterange)
ax2.XTickLabel = '';
ax2.YLabel.String = 'CN';
ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'off';
ax2.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax2.XAxis.TickDirection = 'out';
ax2.YAxis.MinorTick = 'off';
ax2.YAxis.TickDirection = 'out';

title(title_string{2},'FontSize',15)

linkaxes([ax1 ax2],'x')