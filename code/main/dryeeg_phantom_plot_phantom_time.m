%% Phantom Figure TimeSeries
clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m

rmappath=dryeeg_result_path(repoRoot, 'phantom_plot_phantom_time');
if ~exist(rmappath,"file")
    mkdir(rmappath)
end
files={...
    dryeeg_data_path(repoRoot, 'PhantomComp/archive/P007_data_stim_average_240213_1817.mat');...
    dryeeg_data_path(repoRoot, 'PhantomComp/archive/P002_data_stim_average_240213_1759.mat');...
    dryeeg_data_path(repoRoot, 'PhantomComp/archive/Phantom_data_stim_average_240209_1759.mat')};

for bb=1:length(files)
    orig= load(files{bb});
    if bb==1
        Patient=orig.F;
    elseif bb==2
        Patient2=orig.F;
    elseif bb==3
        Phantom=orig.F;
    end

end

EEGabs=abs(Phantom*-10e6);
EMGabs= abs(Patient*10e6);
MEGabs=abs(Patient2*10e6);

cls=mandrill(3);

for pp=1:length(cls)
    c{pp}=cls(pp,:);
end

c{2}=[0.4,0.4,0.4];

time=orig.Time;
h1=figure;

mamp=mean(MEGabs-EEGabs);
stdamp=std(MEGabs);
semamp= std(MEGabs)/ sqrt (size(MEGabs,1));
plot(time,mamp,'Color',c{1},'LineWidth',2)
hold on
ciplot(mamp+semamp,mamp-semamp,time,c{1})
ans.EdgeColor='none';
ans.FaceAlpha=0.1;

mamp=mean(EMGabs);
stdamp=std(EMGabs);
semamp= std(EMGabs)/ sqrt (size(EMGabs,1));
plot(time,mamp,'Color',c{3},'LineWidth',2)
hold on
ciplot(mamp+semamp,mamp-semamp,time,c{3})
ans.EdgeColor='none';
ans.FaceAlpha=0.1;

EEGabs(7,:)=[];

mamp=mean(EEGabs);
stdamp=std(EEGabs);
semamp= std(EEGabs)/ sqrt (size(EEGabs,1));
plot(time,mamp,'Color',c{2},'LineWidth',2)
hold on
ciplot(mamp+semamp,mamp-semamp,time,c{2})
ans.EdgeColor='none';
ans.FaceAlpha=0.1;

xline(0,'--','Linewidth',1.5);

legend('Patient #2','','Patient #7','','Phantom','')

ax=gca;
ax.XGrid='on';
ax.YGrid='on';
ax.XLabel.String='Time (ms)';
ax.XLim=[-0.03 0.15];
ax.YLim=[3 1000];
ax.YScale='log';
ax.YLabel.String='Log Amplitude (ÂµV)';
ax.Legend.Box='off';
ax.TickDir='out';
ax.Box='off';
ax.FontSize=10;
fig=gcf;
fig.Color='white';
ax=gca;
ax.TickDir="out";
ax.Box="off";
ax.XGrid="on";
ax.YGrid="on";
ax.FontSize=14;
ax.FontWeight="bold";

% saveas(h1,[rmappath,'PhantomvsPatients_logscale','.fig']);
% saveas(h1,[rmappath,'PhantomvsPatients_logscale','.png']);


