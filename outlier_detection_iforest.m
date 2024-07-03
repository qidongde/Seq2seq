%% Formatting CN data as 1/second
% Missing values ​​are replaced by NaN
clear;

load('UHSAS_2016.mat');
% Dp_ = mean(Dp_bounds,1)';

% Select one time period
t1 = datetime(2016,1,5,0,0,0);
t2 = datetime(2017,1,1,0,0,0);
t = t1:seconds(1):t2;
uhsas_time_second=t';

time_filter = and(Time_UHSAS>datenum(t1),Time_UHSAS<datenum(t2));
N_selected = N(time_filter);
Time_selected = Time_UHSAS(time_filter);

data_position = seconds(datetime(Time_selected,'ConvertFrom','datenum')-t1)+1;
data_position = round(data_position,TieBreaker="minusinf");
uhsas_cn_second = NaN(numel(t),1);
uhsas_cn_second(data_position,:) = N_selected;

% uhsas_size_dist_second = NaN(numel(t),99);
% uhsas_size_dist_second(data_position,:)=dN_dlogDp;

save( 'uhsas_secondly_selected_time_2016.mat', 'uhsas_cn_second', "uhsas_time_second");

%% step1 iforest data preprocess
clear;

load('uhsas_secondly_selected_time_2016.mat');

CN_1min.ave = movmean(uhsas_cn_second,61,'omitnan');
CN_1min.std = movstd(uhsas_cn_second,61,'omitnan');
CN_1min.median = movmedian(uhsas_cn_second,61,'omitnan');
uhsas_CN_1sContinuous.CN_1min_moving = CN_1min;

CN_2min.ave = movmean(uhsas_cn_second,121,'omitnan');
CN_2min.std = movstd(uhsas_cn_second,121,'omitnan');
CN_2min.median = movmedian(uhsas_cn_second,121,'omitnan');
uhsas_CN_1sContinuous.CN_2min_moving = CN_2min;

CN_5min.ave = movmean(uhsas_cn_second,301,'omitnan');
CN_5min.std = movstd(uhsas_cn_second,301,'omitnan');
CN_5min.median = movmedian(uhsas_cn_second,301,'omitnan');
uhsas_CN_1sContinuous.CN_5min_moving = CN_5min;

CN_30min.ave = movmean(uhsas_cn_second,1801,'omitnan');
CN_30min.std = movstd(uhsas_cn_second,1801,'omitnan');
CN_30min.median = movmedian(uhsas_cn_second,1801,'omitnan');
uhsas_CN_1sContinuous.CN_30min_moving = CN_30min;

save('uhsas_secondly_selected_time_2016.mat','uhsas_CN_1sContinuous','-append')

%% Step1 iforest filter
clear;

load('uhsas_secondly_selected_time_2016.mat');

% Feature engineering:
% CN_raw – CN_30minMedian
deviate = uhsas_cn_second-uhsas_CN_1sContinuous.CN_30min_moving.median;
% CN_1minStd/max(CN_1minAve,300)
SARatio1min = uhsas_CN_1sContinuous.CN_1min_moving.std./max(300,uhsas_CN_1sContinuous.CN_1min_moving.ave); 
% CN_2minStd/max(CN_2minAve,300)
SARatio2min = uhsas_CN_1sContinuous.CN_2min_moving.std./max(300,uhsas_CN_1sContinuous.CN_2min_moving.ave); 
% (CN_raw – CN_1minMedian)/max(CN_1minMedian,300)
mu1min = (uhsas_cn_second-uhsas_CN_1sContinuous.CN_1min_moving.median)./max(300,uhsas_CN_1sContinuous.CN_1min_moving.median); 
% (CN_raw – CN_5minMedian)/max(CN_5minMedian,300)
mu5min = (uhsas_cn_second-uhsas_CN_1sContinuous.CN_5min_moving.median)./max(300,uhsas_CN_1sContinuous.CN_5min_moving.median); 

% Five features
sample = [deviate SARatio1min SARatio2min mu1min mu5min]; 

% remove missing data
ID_valid = find(~isnan(deviate));
sample = sample(ID_valid,:); 

% get the anomaly score from iforest method
[~,~,scores1] = iforest(sample); 
% Hyperparameter. Plot the histogram to find an appropriate threshold.
score_threshold = 0.55; 

iforest_label = NaN(length(uhsas_cn_second),1);
scores = NaN(length(uhsas_cn_second),1);

iforest_label(ID_valid(scores1>=score_threshold)) = 1;
iforest_label(ID_valid(scores1<score_threshold)) = 0;
scores(ID_valid) = scores1;

% 1: polluted identified by Step-1
uhsas_CN_1sContinuous.iForest.AnomalyScore = scores;
uhsas_CN_1sContinuous.iForest.isAnomaly = iforest_label; 

save('uhsas_secondly_selected_time_2016.mat','uhsas_CN_1sContinuous','-append')

%% Score hist --> determine the threshold
clear;

load('uhsas_secondly_selected_time_2016.mat');

y1 = uhsas_CN_1sContinuous.iForest.AnomalyScore;
ID_valid = find(~isnan(y1));
y1 = y1(ID_valid,:); 


x = linspace(0,1,50);
y2 = NaN(size(x));

for i = 1:numel(x)
    y2(i) = sum(y1>x(i))/numel(y1)*100
end

yyaxis left;
plot(x,y2);
ylim([0,100]); 
ylabel('Percentage');

yyaxis right; 
nbins = 50;
histogram(y1,nbins);
% title('Title');
xlabel('Anomaly score');
ylabel('Number'); 
grid on;

%% step1 result visualization
clear;

load('uhsas_secondly_selected_time_2016.mat');
scores = uhsas_CN_1sContinuous.iForest.AnomalyScore;
threshold = 0.55;

t = datenum(uhsas_time_second);
x1 = t(scores<threshold);
y1 = uhsas_cn_second(scores<threshold);
x2 = t(scores>=threshold);
y2 = uhsas_cn_second(scores>=threshold);

fig = figure;
set(fig,'Color','w','Position',[100 100 1000 300]);
scatter(x1,y1,'.');
title('CN');

% xlabel('Date');
ylabel('CN');
datetick('x','yyyy-mm-dd');

hold on;
scatter(x2,y2,'.');
xlim([t(1),t(end)]);
% set(gca,'xtick',t(1):5:t(end));
hold off;

%% Step2 iforest data preprocess
clear;

load('uhsas_secondly_selected_time_2016.mat');

uhsas_CN_1sContinuous_2.cn = uhsas_cn_second;
uhsas_CN_1sContinuous_2.cn(uhsas_CN_1sContinuous.iForest.isAnomaly~=0) = NaN;

CN_10min.ave = movmean(uhsas_CN_1sContinuous_2.cn,601,'omitnan');
CN_10min.std = movstd(uhsas_CN_1sContinuous_2.cn,601,'omitnan');
CN_10min.median = movmedian(uhsas_CN_1sContinuous_2.cn,601,'omitnan');
uhsas_CN_1sContinuous_2.CN_10min_moving = CN_10min;

CN_30min.ave = movmean(uhsas_CN_1sContinuous_2.cn,1801,'omitnan');
CN_30min.std = movstd(uhsas_CN_1sContinuous_2.cn,1801,'omitnan');
CN_30min.median = movmedian(uhsas_CN_1sContinuous_2.cn,1801,'omitnan');
uhsas_CN_1sContinuous_2.CN_30min_moving = CN_30min;

CN_60min.ave = movmean(uhsas_CN_1sContinuous_2.cn,3601,'omitnan');
CN_60min.std = movstd(uhsas_CN_1sContinuous_2.cn,3601,'omitnan');
CN_60min.median = movmedian(uhsas_CN_1sContinuous_2.cn,3601,'omitnan');
uhsas_CN_1sContinuous_2.CN_60min_moving = CN_60min;

save('uhsas_secondly_selected_time_2016.mat','uhsas_CN_1sContinuous_2','-append')

%% Step2 iforest filter

SARatio10min = uhsas_CN_1sContinuous_2.CN_10min_moving.std./max(300,uhsas_CN_1sContinuous_2.CN_10min_moving.ave);
SARatio30min = uhsas_CN_1sContinuous_2.CN_30min_moving.std./max(300,uhsas_CN_1sContinuous_2.CN_30min_moving.ave);
SARatio60min = uhsas_CN_1sContinuous_2.CN_60min_moving.std./max(300,uhsas_CN_1sContinuous_2.CN_60min_moving.ave);

mu10min = (uhsas_CN_1sContinuous_2.cn-uhsas_CN_1sContinuous_2.CN_10min_moving.median)./max(300,uhsas_CN_1sContinuous_2.CN_10min_moving.median);
mu30min = (uhsas_CN_1sContinuous_2.cn-uhsas_CN_1sContinuous_2.CN_30min_moving.median)./max(300,uhsas_CN_1sContinuous_2.CN_30min_moving.median);
mu60min = (uhsas_CN_1sContinuous_2.cn-uhsas_CN_1sContinuous_2.CN_60min_moving.median)./max(300,uhsas_CN_1sContinuous_2.CN_60min_moving.median);

sample = [SARatio10min SARatio30min SARatio60min mu10min mu30min mu60min];

ID_valid = find(~isnan(mu10min));
sample = sample(ID_valid,:);

[~,~,scores2] = iforest(sample);
score_threshold = 0.60; 

tf = NaN(length(uhsas_time_second),1);
scores = NaN(length(uhsas_time_second),1);

tf(ID_valid(scores2>=score_threshold)) = 1;
tf(ID_valid(scores2<score_threshold)) = 0;
scores(ID_valid) = scores2;

uhsas_CN_1sContinuous_2.iForest.isAnomaly = tf;
uhsas_CN_1sContinuous_2.iForest.AnomalyScore = scores;

save('uhsas_secondly_selected_time_2016.mat','uhsas_CN_1sContinuous_2','-append')

%% Score hist --> determine the threshold
clear;

load('uhsas_secondly_selected_time_2016.mat');

y1 = uhsas_CN_1sContinuous_2.iForest.AnomalyScore;
ID_valid = find(~isnan(y1));
y1 = y1(ID_valid,:); 


x = linspace(0,1,50);
y2 = NaN(size(x));

for i = 1:numel(x)
    y2(i) = sum(y1>x(i))/numel(y1)*100
end

yyaxis left;
plot(x,y2);
ylim([0,100]); 
ylabel('Percentage');

yyaxis right; 
nbins = 50;
histogram(y1,nbins);
% title('Title');
xlabel('Anomaly score');
ylabel('Number'); 
grid on;
%% step2 result visualization
clear;

load('uhsas_secondly_selected_time_2016.mat');
scores_1 = uhsas_CN_1sContinuous.iForest.AnomalyScore;
scores_2 = uhsas_CN_1sContinuous_2.iForest.AnomalyScore;
threshold_1 = 0.55;
threshold_2 = 0.6;

t = datenum(uhsas_time_second);
normal_label = and(scores_1<threshold_1,scores_2<threshold_2);
x1 = t(normal_label);
y1 = uhsas_cn_second(normal_label);
filter1 = scores_1>=threshold_1;
x2 = t(filter1);
y2 = uhsas_cn_second(filter1);
filter2 = and(scores_1<threshold_1,scores_2>=threshold_2);
x3 = t(filter2);
y3 = uhsas_cn_second(filter2);

fig = figure;
set(fig,'Color','w','Position',[100 100 1000 400]);
scatter(x1,y1,'.');
title('CN');

% xlabel('Date');
ylabel('CN');
datetick('x','yyyy-mm-dd');

hold on;
scatter(x2,y2,'.');
xlim([t(1),t(end)]);
% set(gca,'xtick',t(1):5:t(end));

hold on;
scatter(x3,y3,'.');
% xlim([t(1),t(end)]);
ylim([0,2000]); 
% set(gca,'xtick',t(1):5:t(end));
legend('Normal','step1','step2');
hold off;

%% Step3 data preprocess
clear;

load('uhsas_secondly_selected_time_2016.mat');
uhsas_CN_1sContinuous_3.cn = uhsas_CN_1sContinuous_2.cn;
uhsas_CN_1sContinuous_3.cn(uhsas_CN_1sContinuous_2.iForest.isAnomaly~=0) = NaN;

CN_60min.ave = movmean(uhsas_CN_1sContinuous_3.cn,3601,'omitnan');
CN_60min.std = movstd(uhsas_CN_1sContinuous_3.cn,3601,'omitnan');
CN_60min.median = movmedian(uhsas_CN_1sContinuous_3.cn,3601,'omitnan');
uhsas_CN_1sContinuous_3.CN_60min_moving = CN_60min;

% Identify the high-CN variation events
std_threshold = 75;
id_high_std = find(uhsas_CN_1sContinuous_3.CN_60min_moving.std > std_threshold);
% Identify the discontinuous ids which represent the borders between different high-variation events
id_diff = diff(id_high_std); 

event_id = [[id_high_std(1);id_high_std(find(id_diff>1)+1)] [id_high_std(id_diff>1);id_high_std(end)]]; % two columns represent the start and end id of each event

event_id_combine = event_id;
% combine adjacent events with time gap less than 1.5h
j = 1;
while j<size(event_id_combine,1) 
    if event_id_combine(j+1,1)-event_id_combine(j,2) < 90*60
        column1 = event_id_combine(:,1);
        column1(j+1) = [];
        column2 = event_id_combine(:,2);
        column2(j) = [];
        event_id_combine = [column1 column2];
    else
        j = j+1;
    end
end

% number of events with high std
event_num = size(event_id_combine,1); 

% and assign pollutionflag
uhsas_CN_1sContinuous_3.spikyflag = zeros(size(uhsas_CN_1sContinuous_3.cn));
for i = 1:event_num
    uhsas_CN_1sContinuous_3.spikyflag(event_id_combine(i,1):event_id_combine(i,2)) = 1;
end
uhsas_CN_1sContinuous_3.spikyflag(isnan(uhsas_CN_1sContinuous_3.cn)) = NaN;

save('uhsas_secondly_selected_time_2016.mat','uhsas_CN_1sContinuous_3','-append')

%% Score hist --> determine the threshold
clear;

load('uhsas_secondly_selected_time.mat');

y1 = uhsas_CN_1sContinuous_3.CN_60min_moving.std;
ID_valid = find(~isnan(y1));
y1 = y1(ID_valid,:); 


x = linspace(0,max(y1),500);
y2 = NaN(size(x));

for i = 1:numel(x)
    y2(i) = sum(y1>x(i))/numel(y1)*100
end

yyaxis left;
plot(x,y2);
ylim([0,100]); 
ylabel('Percentage');

yyaxis right; 
nbins = 50;
histogram(y1,nbins);
% title('Title');
xlim([0,200]);
xlabel('60min moving std');
ylabel('Number'); 
grid on;
%% step3 result visualization
clear;

load('uhsas_secondly_selected_time_2016.mat');
scores_1 = uhsas_CN_1sContinuous.iForest.AnomalyScore;
scores_2 = uhsas_CN_1sContinuous_2.iForest.AnomalyScore;
threshold_1 = 0.55;
threshold_2 = 0.6;

t = datenum(uhsas_time_second);
spike_filter = (uhsas_CN_1sContinuous_3.spikyflag==1);
normal_label = and(scores_1<threshold_1,scores_2<threshold_2) & ~spike_filter;
x1 = t(normal_label);
y1 = uhsas_cn_second(normal_label);
filter1 = and(scores_1>=threshold_1,~spike_filter);
x2 = t(filter1);
y2 = uhsas_cn_second(filter1);
filter2 = and(scores_1<threshold_1,scores_2>=threshold_2) & ~spike_filter;
x3 = t(filter2);
y3 = uhsas_cn_second(filter2);
filter3 = and(scores_1<threshold_1,scores_2<threshold_2) & spike_filter;
x4 = t(filter3);
y4 = uhsas_cn_second(filter3);

fig = figure;
set(fig,'Color','w','Position',[100 100 1000 400]);
scatter(x1,y1,'.');
title('CN');

% xlabel('Date');
ylabel('CN');
datetick('x','yyyy-mm-dd');

hold on;
scatter(x2,y2,'.');

hold on;
scatter(x3,y3,'.');

hold on;
scatter(x4,y4,'.');
xlim([t(1),t(end)]);
ylim([0,2000]); 
% set(gca,'xtick',t(1):5:t(end));
legend('Normal','step1','step2','step3');
hold off;

%% clean hourly data
clear;

load('uhsas_secondly_selected_time_2016.mat');
load('UHSAS_2016.mat');
Dp_ = mean(Dp_bounds,1)';

% uhsas size distribution data secondly
t1 = datetime(2016,1,5,0,0,0);
t2 = datetime(2017,1,1,0,0,0);
time_filter = and(Time_UHSAS>datenum(t1),Time_UHSAS<datenum(t2));
dN_dlogDp_selected = dN_dlogDp(time_filter,:);
Time_selected = Time_UHSAS(time_filter);

data_position = seconds(datetime(Time_selected,'ConvertFrom','datenum')-t1)+1;
data_position = round(data_position,TieBreaker="minusinf");
uhsas_sd_second = NaN(numel(uhsas_time_second),99);
uhsas_sd_second(data_position,:) = dN_dlogDp_selected;

% Filter
ind = find(uhsas_CN_1sContinuous_3.spikyflag == 0);
uhsas_sd_second_clean = uhsas_sd_second(ind,:);
uhsas_cn_second_clean = uhsas_cn_second(ind,:);
uhsas_time_second_clean = uhsas_time_second(ind,:);

% Hourly median calculate
time_vec = datevec(uhsas_time_second_clean);
[G,year,month,day,hour] = findgroups(time_vec(:,1),time_vec(:,2),time_vec(:,3),time_vec(:,4));
time_vec_hourly = [year month day hour];
uhsas_timevec_hour = [time_vec_hourly,zeros(size(time_vec_hourly,1),2)];

% Hourly total CN median
uhsas_cn_hour = splitapply(@median,uhsas_cn_second_clean,G);

% Hourly size distribution median
uhsas_sd_hour = [];
for i=1:99
    tmp = splitapply(@median,uhsas_sd_second_clean(:,i),G);
    uhsas_sd_hour = [uhsas_sd_hour,tmp];
end

% Data organization
t1 = datetime(2016,1,5,0,0,0);
t2 = datetime(2017,1,1,0,0,0);
t = t1:hours(1):t2;
uhsas_time_hour=t';

data_position = hours(datetime(uhsas_timevec_hour)-t1)+1;
data_position = round(data_position,TieBreaker="minusinf");
uhsas_sd_hour_all = NaN(numel(uhsas_time_hour),99);
uhsas_sd_hour_all(data_position,:) = uhsas_sd_hour;

uhsas_cn_hour_all = NaN(numel(uhsas_time_hour),1);
uhsas_cn_hour_all(data_position,:) = uhsas_cn_hour;

save('uhsas_clean_data_hourly_2016.mat','uhsas_sd_hour_all', 'uhsas_cn_hour_all','uhsas_time_hour', 'Dp_');
%% visualization of clean data
clear;

load('uhsas_clean_data_hourly_2016.mat');

daterange = [datenum(2016,1,1,0,0,0) datenum(2017,1,1,0,0,0)];


title_string = {'Size distribution','CN'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1800 700])


% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.3]);

time_uhsas_heated = datenum(uhsas_time_hour)
idx_uhsas_heated = find(time_uhsas_heated>=daterange(1) & time_uhsas_heated<daterange(2));

x_all = time_uhsas_heated(idx_uhsas_heated);
y1 = uhsas_sd_hour_all(idx_uhsas_heated,:);

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

y2 = uhsas_cn_hour_all(idx_uhsas_heated,:);
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
