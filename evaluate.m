%% train test split
clear;

load('Dataset_input_output_2014_2023.mat');

test_filter = find(time_selected_vector(:,1)==2022);
train_filter = find(time_selected_vector(:,1)~=2022);

test_input = input_selected(test_filter,:);
test_output = output_selected(test_filter,:);
test_time = time_selected(test_filter,:);
test_time_vec = time_selected_vector(test_filter,:);

train_input = input_selected(train_filter,:);
train_output = output_selected(train_filter,:);
train_time = time_selected(train_filter,:);
train_time_vec = time_selected_vector(train_filter,:);

save('Split_Dataset.mat','train_input','train_output','train_time', ...
    "train_time_vec",'test_input','test_output','test_time', ...
    'test_time_vec','Dp_');

%% visualization
figure('Position',[100 100 1400 700])
for i = 1:30
    subplot(5,6,i)
    plot(10.^train_y_true_final(i+30000,:));
    hold on
    plot(10.^train_y_pre_final(i+30000,:));
    ylim([0 1000])
end

%%
figure;
for i=1:30
    subplot(5,6,i)
    scatter(train_y_pre_final(1,:),train_y_pre_final(100*i,:))
end
%% test
figure('Position',[100 100 1400 700])
for i = 1:30
    subplot(5,6,i)
    plot(10.^test_y_true_final(i+3100,:));
    hold on
    plot(10.^test_y_pre_final(i+3100,:));
    ylim([0 1000])
end

%%
figure;
for i=1:30
    subplot(5,6,i)
    scatter(test_y_pre_final(1,:),test_y_pre_final(100*i,:))
end

%% train compare
clear;

load('Split_Dataset.mat');
load('LSTM_FC_50_Result_.mat');

t1 = datetime(2014,1,1,0,0,0);
t2 = datetime(2024,1,1,0,0,0);
t = t1:hours(1):t2;
hourly_time=t';

train_label = hours(datetime(train_time_vec)-t1)+1;
train_true = NaN(numel(t),99);
train_predict = NaN(numel(t),99);

train_true(train_label,:)=10.^train_y_true_final-1;
train_predict(train_label,:)=10.^train_y_pre_final-1;

test_label = hours(datetime(test_time_vec)-t1)+1;
test_true = NaN(numel(t),99);
test_predict = NaN(numel(t),99);

test_true(test_label,:)=10.^test_y_true_final-1;
test_predict(test_label,:)=10.^test_y_pre_final-1;


train_daterange = [datenum(2014,1,1,0,0,0) datenum(2024,1,1,0,0,0)];

title_string = {'Train Ground Truth','Train Prediction'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1200 600])


% GT heatmap
ax1 = axes('Position',[0.08 0.08 0.80 0.25]);

time_stamp = datenum(hourly_time);
train_time_label = find(time_stamp>=train_daterange(1) & ...
    time_stamp<train_daterange(2));

train_time = hourly_time(train_time_label);
train_true_heated = train_true(train_time_label,:);
PC1 = pcolor(train_time,Dp_,train_true_heated');

set(PC1,'EdgeColor','none')
caxis([0 1200])
ax1.YScale = 'log';
set(ax1,'FontSize',12)
% ax1.Xlim(train_daterange)
datetick('x','yyyy-mm-dd')

ax1.YLabel.String = 'D_p (nm)';
ax1.YLabel.FontSize = 15;
ax1.XAxis.MinorTick = 'off';
% ax1.XAxis.MinorTickValues = [train_daterange(1):30:train_daterange(2)];
ax1.XAxis.TickDirection = 'out';
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.TickDirection = 'out';
ax1.YLim = [60 1000];
title(title_string{1},'FontSize',15)

h1 = colorbar('v');
h1.Limits = [0 1200];
colormap(jet)
h1.Position = [0.89 0.08 0.015 0.25];
h1.FontSize = 12;
h1.Label.String = 'd{\itN}/d{\itlogD_p} (cm^{–3})';
h1.Label.FontSize = 15;
h1.Ticks = [0:300:1200];
h1.TickDirection = 'out';
h1.TickLength = 0.03;

% Prediction heatmap
ax2 = axes('Position',[0.08 0.48 0.80 0.25]);

train_pre_heated = train_predict(train_time_label,:);
PC2 = pcolor(train_time,Dp_,train_pre_heated');

set(PC2,'EdgeColor','none')
caxis([0 1200])
ax2.YScale = 'log';
set(ax2,'FontSize',12)
% ax1.Xlim(test_daterange)
datetick('x','yyyy-mm-dd')

% ax2.YLabel.String = 'D_p (nm)';
ax2.XTickLabel = '';
ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'off';
% ax1.XAxis.MinorTickValues = [train_daterange(1):30:train_daterange(2)];
ax2.XAxis.TickDirection = 'out';
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.TickDirection = 'out';
ax2.YLim = [60 1000];
title(title_string{2},'FontSize',15)

h2 = colorbar('v');
h2.Limits = [0 1200];
colormap(jet)
h2.Position = [0.89 0.48 0.015 0.25];
h2.FontSize = 12;
% h2.Label.String = 'd{\itN}/d{\itlogD_p} (cm^{–3})';
h2.Label.FontSize = 15;
h2.Ticks = [0:300:1200];
h2.TickDirection = 'out';
h2.TickLength = 0.03;

linkaxes([ax1 ax2],'x')


%% test compare
clear;

load('Split_Dataset.mat');
load('Seq2Seq_200_Result_.mat');

t1 = datetime(2014,1,1,0,0,0);
t2 = datetime(2024,1,1,0,0,0);
t = t1:hours(1):t2;
hourly_time=t';

train_label = hours(datetime(train_time_vec)-t1)+1;
train_true = NaN(numel(t),99);
train_predict = NaN(numel(t),99);

train_true(train_label,:)=10.^train_y_true_final-1;
train_predict(train_label,:)=10.^train_y_pre_final-1;

test_label = hours(datetime(test_time_vec)-t1)+1;
test_true = NaN(numel(t),99);
test_predict = NaN(numel(t),99);

test_true(test_label,:)=10.^test_y_true_final-1;
test_predict(test_label,:)=10.^test_y_pre_final-1;


test_daterange = [datenum(2022,1,1,0,0,0) datenum(2023,1,1,0,0,0)];

title_string = {'Test Ground Truth','Test Prediction'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1400 600])


% GT heatmap
ax1 = axes('Position',[0.08 0.08 0.80 0.25]);

time_stamp = datenum(hourly_time);
test_time_label = find(time_stamp>=test_daterange(1) & ...
    time_stamp<test_daterange(2));

test_time = hourly_time(test_time_label);
test_true_heated = test_true(test_time_label,:);
PC1 = pcolor(test_time,Dp_,test_true_heated');

set(PC1,'EdgeColor','none')
caxis([0 1200])
ax1.YScale = 'log';
set(ax1,'FontSize',12)
% ax1.Xlim(train_daterange)
datetick('x','yyyy-mm-dd')

ax1.YLabel.String = 'D_p (nm)';
ax1.YLabel.FontSize = 15;
ax1.XAxis.MinorTick = 'off';
% ax1.XAxis.MinorTickValues = [train_daterange(1):30:train_daterange(2)];
ax1.XAxis.TickDirection = 'out';
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.TickDirection = 'out';
ax1.YLim = [60 1000];
title(title_string{1},'FontSize',15)

h1 = colorbar('v');
h1.Limits = [0 1200];
colormap(jet)
h1.Position = [0.89 0.08 0.015 0.25];
h1.FontSize = 12;
h1.Label.String = 'd{\itN}/d{\itlogD_p} (cm^{–3})';
h1.Label.FontSize = 15;
h1.Ticks = [0:300:1200];
h1.TickDirection = 'out';
h1.TickLength = 0.03;

% Prediction heatmap
ax2 = axes('Position',[0.08 0.48 0.80 0.25]);

test_pre_heated = test_predict(test_time_label,:);
PC2 = pcolor(test_time,Dp_,test_pre_heated');

set(PC2,'EdgeColor','none')
caxis([0 1200])
ax2.YScale = 'log';
set(ax2,'FontSize',12)
% ax1.Xlim(test_daterange)
datetick('x','yyyy-mm-dd')

% ax2.YLabel.String = 'D_p (nm)';
ax2.XTickLabel = '';
ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'off';
% ax1.XAxis.MinorTickValues = [train_daterange(1):30:train_daterange(2)];
ax2.XAxis.TickDirection = 'out';
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.TickDirection = 'out';
ax2.YLim = [60 1000];
title(title_string{2},'FontSize',15)

h2 = colorbar('v');
h2.Limits = [0 1200];
colormap(jet)
h2.Position = [0.89 0.48 0.015 0.25];
h2.FontSize = 12;
% h2.Label.String = 'd{\itN}/d{\itlogD_p} (cm^{–3})';
h2.Label.FontSize = 15;
h2.Ticks = [0:300:1200];
h2.TickDirection = 'out';
h2.TickLength = 0.03;

linkaxes([ax1 ax2],'x')

%% train compare
clear;

load('Split_Dataset.mat');
load('LSTM_FC_50_Result_.mat');

train_true=10.^train_y_true_final-1;
train_predict=10.^train_y_pre_final-1;

y1 = mean(train_true,1);
y2 = mean(train_predict,1);

load('LSTM_FC_100_Result_.mat');
train_predict=10.^train_y_pre_final-1;
y3 = mean(train_predict,1);

load('LSTM_FC_200_Result_.mat');
train_predict=10.^train_y_pre_final-1;
y4 = mean(train_predict,1);

load('Seq2Seq_50_Result_.mat');
train_predict=10.^train_y_pre_final-1;
y5 = mean(train_predict,1);

semilogx(Dp_,y1,Dp_,y2,Dp_,y3,Dp_,y4,Dp_,y5);
legend('Train Ground Truth Average','50 Train Prediction Average', ...
    '100 Train Prediction Average','200 Train Prediction Average', ...
    '50 Train Prediction Average Seq2Seq');
% title('T-x');
xlabel('Dp');
ylabel('d{\itN}/d{\itlogD_p} (cm^{–3})');

%% test compare
clear;

load('Split_Dataset.mat');
load('Seq2Seq_50_Result_.mat');

test_true=10.^test_y_true_final-1;
test_predict=10.^test_y_pre_final-1;

y1 = mean(test_true,1);
y2 = mean(test_predict,1);

load('Seq2Seq_100_Result_.mat');
test_predict=10.^test_y_pre_final-1;
y3 = mean(test_predict,1);

load('Seq2Seq_200_Result_.mat');
test_predict=10.^test_y_pre_final-1;
y4 = mean(test_predict,1);

load('Seq2Seq_50_Result_.mat');
test_predict=10.^test_y_pre_final-1;
y5 = mean(test_predict,1);

semilogx(Dp_,y1,Dp_,y2,Dp_,y3,Dp_,y4);
legend('Test Ground Truth Average','50 Test Prediction Average', ...
    '100 Test Prediction Average','200 Test Prediction Average');
% title('T-x');
xlabel('Dp');
ylabel('d{\itN}/d{\itlogD_p} (cm^{–3})');

%% y=x compare
clear;

load('Split_Dataset.mat');
load('LSTM_FC_200_Result_.mat');

test_true=10.^test_y_true_final-1;
test_predict=10.^test_y_pre_final-1;

x = reshape(test_true,1,[]);
y = reshape(test_predict,1,[]);

% R-squared
SSres = sum((y - x).^2); 
SStot = sum((x - mean(x)).^2); 
R_squared = 1 - (SSres / SStot); 

% cosine value
cos_seq = NaN(numel(test_time),1);
for i=1:numel(test_time)
    cos_seq(i,:)=dot(test_true(i,:),test_predict(i,:))/(norm(test_true(i,:))*norm(test_predict(i,:)));
end
cos_avg = mean(cos_seq);

% rmse
RMSE = sqrt(mean((y-x).^2));

% 60-300
x1 = reshape(test_true(:,1:58),1,[]);
y1 = reshape(test_predict(:,1:58),1,[]);

% R-squared
SSres1 = sum((y1 - x1).^2); 
SStot1 = sum((x1 - mean(x1)).^2); 
R_squared1 = 1 - (SSres1 / SStot1); 

% cosine value
cos_seq1 = NaN(numel(test_time),1);
for i=1:numel(test_time)
    cos_seq1(i,:)=dot(test_true(i,1:58),test_predict(i,1:58))/(norm(test_true(i,1:58))*norm(test_predict(i,1:58)));
end
cos_avg1 = mean(cos_seq1);

% rmse
RMSE1 = sqrt(mean((y1-x1).^2));

% 300-1000
x2 = reshape(test_true(:,59:99),1,[]);
y2 = reshape(test_predict(:,59:99),1,[]);

% R-squared
SSres2 = sum((y2 - x2).^2); 
SStot2 = sum((x2 - mean(x2)).^2); 
R_squared2 = 1 - (SSres2 / SStot2); 

% cosine value
cos_seq2 = NaN(numel(test_time),1);
for i=1:numel(test_time)
    cos_seq2(i,:)=dot(test_true(i,59:99),test_predict(i,59:99))/(norm(test_true(i,59:99))*norm(test_predict(i,59:99)));
end
cos_avg2 = mean(cos_seq2,'omitnan');

% rmse
RMSE2 = sqrt(mean((y2-x2).^2));


set(gcf,'Units','centimeters','Position',[6 6 14 13]);
scatter(x,y,'.');
set(gca,'Xlim',[0,1500],'Ylim',[0,1500],'XTick',[0:300:1500], ...
    'YTick',[0:300:1500]);
xlabel('Test Ground Truth'); 
ylabel('Test Prediction'); 

hold on;
h1=refline(1,0);
set(h1,'color','black','linewidth',1.5);
text(100,1300,'R^2=0.5656');
text(100,1200,'N=553410');


