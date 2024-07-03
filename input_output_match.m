%% Data correction
clear;

load('uhsas_clean_data_hourly_2016.mat');
% uhsas_sd_hour_all(:,52)=(uhsas_sd_hour_all(:,51)+uhsas_sd_hour_all(:,53))/2;
% uhsas_sd_hour_all(:,27)=(uhsas_sd_hour_all(:,28)+uhsas_sd_hour_all(:,26))/2;

semilogx(Dp_,uhsas_sd_hour_all(6,:),Dp_,uhsas_sd_hour_all(16,:), ...
    Dp_,uhsas_sd_hour_all(26,:),Dp_,uhsas_sd_hour_all(36,:),Dp_, ...
    uhsas_sd_hour_all(46,:), Dp_,uhsas_sd_hour_all(56,:));

legend('Time 1','Time 2','Time 3','Time 4', 'Time 5', 'Time 6');
xlabel('D_p (nm)');
ylabel('d{\itN}/d{\itlogD_p} (cm^{–3})');

%% Data match
clear;

load('uhsas_clean_data_hourly_2016.mat');
load('ENA_ML_input_features_2013_2023_v4.mat');

% data correction
uhsas_sd_hour_all(:,52)=(uhsas_sd_hour_all(:,51)+uhsas_sd_hour_all(:,53))/2;
uhsas_sd_hour_all(:,27)=(uhsas_sd_hour_all(:,28)+uhsas_sd_hour_all(:,26))/2;

t1 = uhsas_time_hour(1,:);
t2 = uhsas_time_hour(end,:);

Time_input = time_traj(:,7);
time_filter = and(Time_input >= datenum(t1),Time_input <= datenum(t2));
Time_input_selected = Time_input(time_filter);
X_raw_selected = X_raw(time_filter);

data_position = hours(datetime(Time_input_selected,'ConvertFrom','datenum')-t1)+1;
data_position = round(data_position,TieBreaker="minusinf");
uhsas_sd_hour_selected = uhsas_sd_hour_all(data_position,:);

% remove missing value
nan_filter = ~any(isnan(uhsas_sd_hour_selected),2);
input_selected = X_raw_selected(nan_filter);
output_selected = uhsas_sd_hour_selected(nan_filter,:);
time_selected = Time_input_selected(nan_filter);

time_selected_vector = datevec(time_selected);
save( 'Dataset_selected_2016.mat', 'input_selected', "output_selected",'time_selected_vector','Dp_');

%% Data visialization
clear;

load('Dataset_selected_2016.mat');

% fill missing value with NaN
t1 = datetime(2016,1,1,0,0,0);
t2 = datetime(2017,1,1,0,0,0);
t = t1:hours(1):t2;
time_hour=t';

time_filter = and(time_selected>=datenum(t1),time_selected<=datenum(t2));
output_selected = output_selected(time_filter,:);
time_selected = time_selected(time_filter);

data_position = hours(datetime(time_selected,'ConvertFrom','datenum')-t1)+1;
data_position = round(data_position,TieBreaker="minusinf");
uhsas_sd_hour = NaN(numel(t),99);
uhsas_sd_hour(data_position,:) = output_selected;

% visualization
daterange = [datenum(2016,1,1,0,0,0) datenum(2017,1,1,0,0,0)];


title_string = {'Size distribution','CN'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1800 700])


% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.3]);

time_uhsas_heated = datenum(time_hour)
idx_uhsas_heated = find(time_uhsas_heated>=daterange(1) & time_uhsas_heated<=daterange(2));

x_all = time_uhsas_heated(idx_uhsas_heated);
y1 = uhsas_sd_hour(idx_uhsas_heated,:);

% Data correction
% y1(:,52)=(y1(:,51)+y1(:,53))/2;
% y1(:,27)=(y1(:,28)+y1(:,26))/2;

PC1 = pcolor(x_all,Dp_,y1');

set(PC1,'EdgeColor','none')
caxis([0 1200])
ax1.YScale = 'log';
set(ax1,'FontSize',12);
xlim(ax1, daterange);
datetick('x','yyyy-mm');
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
