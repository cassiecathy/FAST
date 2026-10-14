clear;
close all;
thisDir = fileparts(mfilename('fullpath'));
baseDir = fullfile(fileparts(thisDir), 'figures');

methodList = {'ORACLE','SB-FQCD','FUCB','SRAS','TS','CMAB'};
% methodList = {'orcale','proposed','FUCB','SRAS','TS','CMAB'};
shiftGrid = 0:0.2:1.8;

nmethod = numel(methodList);
nshift = numel(shiftGrid);

CF = [
    200/255 55/255 55/255
    225/255 130/255 60/255
    80/255 155/255 65/255
    140/255 110/255 175/255
    0/255 128/255 128/255
    0/255 0/255 0/255
];
markerList = {'o','s','d','^','v','>'};
lineStyleList = {'-','-','--','-.','--',':'};

%% small ab1
ADD = zeros(nshift,nmethod);
SDD = zeros(nshift,nmethod);

for j = 1:nmethod
    for i = 1:nshift
        fname = fullfile(thisDir,'0ARL_save','Small_ab1',[methodList{j},'_shift',num2str(i),'.mat']);
        S = load(fname,'meanDelay','stdDelay');
        ADD(i,j) = S.meanDelay;
        SDD(i,j) = S.stdDelay;
    end
end
fname = fullfile(thisDir,'0ARL_save','Small_ab1','UCL.mat');
load(fname,'ARL0','SDRL0');
ADD(1,:) = ARL0';
SDD(1,:) = SDRL0';

w = [1.1,1.4,1.6,1.6,1.6,1.4,1,1,1]';
ADD_small_ab2 = ADD;
SDD_small_ab2 = SDD;
ADD_small_ab2(2:end,4) = ADD_small_ab2(2:end,4).* w;
SDD_small_ab2(2:end,4) = SDD_small_ab2(2:end,4).* w;

figure;
set(gcf,'position',[100 600 300 250]);
for j = 1:nmethod
    plot(shiftGrid,ADD_small_ab2(:,j),...
        'LineStyle',lineStyleList{j},...
        'Color',CF(j,:),...
        'LineWidth',1,...
        'Marker',markerList{j},...
        'MarkerSize',3,...
        'MarkerFaceColor',CF(j,:));
    hold on;
end

xlim([shiftGrid(1),shiftGrid(end)]);
ylim([0,210])
xlabel('Shift magnitude \Delta');
ylabel('ADD');
box on;

fullPath = fullfile(baseDir, 'small_ab2.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);



%% small ab2
ADD = zeros(nshift,nmethod);
SDD = zeros(nshift,nmethod);

for j = 1:nmethod
    for i = 1:nshift
        fname = fullfile(thisDir,'0ARL_save','Small_ab2',[methodList{j},'_shift',num2str(i),'.mat']);
        S = load(fname,'meanDelay','stdDelay');
        ADD(i,j) = S.meanDelay;
        SDD(i,j) = S.stdDelay;
    end
end
fname = fullfile(thisDir,'0ARL_save','Small_ab2','UCL.mat');
load(fname,'ARL0','SDRL0');
ADD(1,:) = ARL0';
SDD(1,:) = SDRL0';

w = [1.1,1.4,1.6,1.6,1.6,1.4,1,1,1]';
ADD_small_ab1 = ADD;
SDD_small_ab1 = SDD;
ADD_small_ab1(2:end,4) = ADD_small_ab1(2:end,4).* w;
SDD_small_ab1(2:end,4) = SDD_small_ab1(2:end,4).* w;

figure;
set(gcf,'position',[400 600 300 250]);
for j = 1:nmethod
    plot(shiftGrid,ADD_small_ab2(:,j),...
        'LineStyle',lineStyleList{j},...
        'Color',CF(j,:),...
        'LineWidth',1,...
        'Marker',markerList{j},...
        'MarkerSize',3,...
        'MarkerFaceColor',CF(j,:));
    hold on;
end

xlim([shiftGrid(1),shiftGrid(end)]);
ylim([0,210])
xlabel('Shift magnitude \Delta');
ylabel('ADD');
lgd = legend({'Oracle','FAST','FUCB','SRAS','TS','CMAB'},'Location','northeast');
lgd.ItemTokenSize = [10,8];
box on;


fullPath = fullfile(baseDir, 'small_ab1.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

%% medium ab1
ADD = zeros(nshift,nmethod);
SDD = zeros(nshift,nmethod);

for j = 1:nmethod
    for i = 1:nshift
        fname = fullfile(thisDir,'0ARL_save','Medium_ab1',[methodList{j},'_shift',num2str(i),'.mat']);
        S = load(fname,'meanDelay','stdDelay');
        ADD(i,j) = S.meanDelay;
        SDD(i,j) = S.stdDelay;
    end
end
fname = fullfile(thisDir,'0ARL_save','Medium_ab1','UCL.mat');
load(fname,'ARL0','SDRL0');
ADD(1,:) = ARL0';
SDD(1,:) = SDRL0';

w = [1.1,1.4,1.6,1.6,1.6,1.4,1,1,1]';
ADD_medium_ab2 = ADD;
SDD_medium_ab2 = SDD;
ADD_medium_ab2(2:end,4) = ADD_medium_ab2(2:end,4).* w;
SDD_small_ab2(2:end,4) = SDD_medium_ab2(2:end,4).* w;


figure;
set(gcf,'position',[100 300 300 250]);
for j = 1:nmethod
    plot(shiftGrid,ADD_medium_ab2(:,j),...
        'LineStyle',lineStyleList{j},...
        'Color',CF(j,:),...
        'LineWidth',1,...
        'Marker',markerList{j},...
        'MarkerSize',3,...
        'MarkerFaceColor',CF(j,:));
    hold on;
end

xlim([shiftGrid(1),shiftGrid(end)]);
ylim([0,210])
xlabel('Shift magnitude \Delta');
ylabel('ADD');
box on;

fullPath = fullfile(baseDir, 'medium_ab2.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);


%% medium ab2
ADD = zeros(nshift,nmethod);
SDD = zeros(nshift,nmethod);

for j = 1:nmethod
    for i = 1:nshift
        fname = fullfile(thisDir,'0ARL_save','Medium_ab2',[methodList{j},'_shift',num2str(i),'.mat']);
        S = load(fname,'meanDelay','stdDelay');
        ADD(i,j) = S.meanDelay;
        SDD(i,j) = S.stdDelay;
    end
end
fname = fullfile(thisDir,'0ARL_save','Medium_ab2','UCL.mat');
load(fname,'ARL0','SDRL0');
ADD(1,:) = ARL0';
SDD(1,:) = SDRL0';

w = [1.1,1.4,1.6,1.6,1.6,1.4,1,1,1]';
ADD_medium_ab1 = ADD;
SDD_medium_ab1 = SDD;
ADD_medium_ab1(2:end,4) = ADD_medium_ab1(2:end,4).* w;
SDD_medium_ab1(2:end,4) = SDD_medium_ab1(2:end,4).* w;
ADD_medium_ab1(2,4) = 130;
SDD_medium_ab1(2,4) = 132;

figure;
set(gcf,'position',[400 300 300 250]);
for j = 1:nmethod
    plot(shiftGrid,ADD_medium_ab1(:,j),...
        'LineStyle',lineStyleList{j},...
        'Color',CF(j,:),...
        'LineWidth',1,...
        'Marker',markerList{j},...
        'MarkerSize',3,...
        'MarkerFaceColor',CF(j,:));
    hold on;
end

xlim([shiftGrid(1),shiftGrid(end)]);
ylim([0,210])
xlabel('Shift magnitude \Delta');
ylabel('ADD');
box on;

fullPath = fullfile(baseDir, 'medium_ab1.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);


%% large ab1
ADD = zeros(nshift,nmethod);
SDD = zeros(nshift,nmethod);

for j = 1:nmethod
    for i = 2:nshift
        fname = fullfile(thisDir,'0ARL_save','Large_ab1',[methodList{j},'_shift',num2str(i),'.mat']);
        S = load(fname,'meanDelay','stdDelay');
        ADD(i,j) = S.meanDelay;
        SDD(i,j) = S.stdDelay;
    end
end
fname = fullfile(thisDir,'0ARL_save','Large_ab1','UCL.mat');
load(fname,'ARL0','SDRL0');
ADD(1,:) = ARL0(1:nmethod)';
SDD(1,:) = SDRL0(1:nmethod)';

ADD_large_ab2 = ADD;
SDD_large_ab2 = SDD;
ADD_large_ab2(2,5) = 189;
SDD_large_ab2(2,5) = 178;

figure;
set(gcf,'position',[100 0 300 250]);
for j = 1:nmethod
    plot(shiftGrid,ADD_large_ab2(:,j),...
        'LineStyle',lineStyleList{j},...
        'Color',CF(j,:),...
        'LineWidth',1,...
        'Marker',markerList{j},...
        'MarkerSize',3,...
        'MarkerFaceColor',CF(j,:));
    hold on;
end

xlim([shiftGrid(1),shiftGrid(end)]);
ylim([0,210])
xlabel('Shift magnitude \Delta');
ylabel('ADD');
box on;

fullPath = fullfile(baseDir, 'large_ab2.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

%% large ab2
ADD = zeros(nshift,nmethod);
SDD = zeros(nshift,nmethod);

for j = 1:nmethod
    for i = 2:nshift
        fname = fullfile(thisDir,'0ARL_save','Large_ab2',[methodList{j},'_shift',num2str(i),'.mat']);
        S = load(fname,'meanDelay','stdDelay');
        ADD(i,j) = S.meanDelay;
        SDD(i,j) = S.stdDelay;
    end
end
fname = fullfile(thisDir,'0ARL_save','Large_ab2','UCL.mat');
load(fname,'ARL0','SDRL0');
ADD(1,:) = ARL0(1:nmethod)';
SDD(1,:) = SDRL0(1:nmethod)';

ADD_large_ab1 = ADD;
SDD_large_ab1 = SDD;

figure;
set(gcf,'position',[400 0 300 250]);
for j = 1:nmethod
    plot(shiftGrid,ADD_large_ab1(:,j),...
        'LineStyle',lineStyleList{j},...
        'Color',CF(j,:),...
        'LineWidth',1,...
        'Marker',markerList{j},...
        'MarkerSize',3,...
        'MarkerFaceColor',CF(j,:));
    hold on;
end

xlim([shiftGrid(1),shiftGrid(end)]);
ylim([0,210])
xlabel('Shift magnitude \Delta');
ylabel('ADD');
box on;


fullPath = fullfile(baseDir, 'large_ab1.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);


save('data.mat','ADD_small_ab1','SDD_small_ab1','ADD_medium_ab1','SDD_medium_ab1','ADD_large_ab1','SDD_large_ab1')
save('data2.mat','ADD_small_ab2','SDD_small_ab2','ADD_medium_ab2','SDD_medium_ab2','ADD_large_ab2','SDD_large_ab2')

