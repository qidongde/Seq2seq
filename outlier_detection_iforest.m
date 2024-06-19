%% Formatting CN data as 1/second
% Missing values ​​are replaced by NaN
clear;

load('UHSAS_2014.mat');
% Dp_ = mean(Dp_bounds,1)';

% Select one time period
t1 = datetime(2014,2,20,0,0,0);
t2 = datetime(2014,3,6,0,0,0);
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

save( 'uhsas_secondly_selected_time.mat', 'uhsas_cn_second', "uhsas_time_second");

%% step1 iforest data preprocess
clear;

load('uhsas_secondly_selected_time.mat');

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

save('uhsas_secondly_selected_time.mat','uhsas_CN_1sContinuous','-append')

%% Step1 iforest filter
clear;

load('uhsas_secondly_selected_time.mat');

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

save('uhsas_secondly_selected_time.mat','uhsas_CN_1sContinuous','-append')

%% Score hist --> determine the threshold
clear;

load('uhsas_secondly_selected_time.mat');

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

load('uhsas_secondly_selected_time.mat');
scores = uhsas_CN_1sContinuous.iForest.AnomalyScore;
threshold = 0.55;

t = datenum(uhsas_time_second);
x1 = t(scores<threshold);
y1 = uhsas_cn_second(scores<threshold);
x2 = t(scores>=threshold);
y2 = uhsas_cn_second(scores>=threshold);

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
set(gca,'xtick',t(1):5:t(end));
hold off;

%% Step2 iforest data preprocess
clear;

load('uhsas_secondly_selected_time.mat');

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

save('uhsas_secondly_selected_time.mat','uhsas_CN_1sContinuous_2','-append')

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

save('uhsas_secondly_selected_time.mat','uhsas_CN_1sContinuous_2','-append')

%% Score hist --> determine the threshold
clear;

load('uhsas_secondly_selected_time.mat');

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

load('uhsas_secondly_selected_time.mat');
scores_1 = uhsas_CN_1sContinuous.iForest.AnomalyScore;
scores_2 = uhsas_CN_1sContinuous_2.iForest.AnomalyScore;
threshold_1 = 0.6;
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
set(gca,'xtick',t(1):5:t(end));

hold on;
scatter(x3,y3,'.');
xlim([t(1),t(end)]);
set(gca,'xtick',t(1):5:t(end));
legend('Normal','step1','step2');
hold off;

%% Step3 data preprocess
clear;

load('uhsas_secondly_selected_time.mat');
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

save('uhsas_secondly_selected_time.mat','uhsas_CN_1sContinuous_3','-append')

%% step3 result visualization
clear;

load('uhsas_secondly_selected_time.mat');
scores_1 = uhsas_CN_1sContinuous.iForest.AnomalyScore;
scores_2 = uhsas_CN_1sContinuous_2.iForest.AnomalyScore;
threshold_1 = 0.6;
threshold_2 = 0.7;

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
set(gca,'xtick',t(1):5:t(end));
legend('Normal','step1','step2','step3');
hold off;