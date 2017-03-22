function ConvertToMat()

global ExportDataSpreadsheet handletext F_samp
[ChannelNumbers,Filenames] = xlsread(ExportDataSpreadsheet,1,'E4:S5000');
F_samp = xlsread(ExportDataSpreadsheet,2,'C3');
[~,ExportDataDirectory,~] = xlsread(ExportDataSpreadsheet,2,'C4');
ConvertMatFlag = ChannelNumbers(:,1);
ChannelNumbers(:,1)=[];
for n=1:size(Filenames,1)
    errorlist=[];
    try
        if ConvertMatFlag(n)
            [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = ImportOneToMat(n,Filenames,ChannelNumbers);
            fname=[Filenames{n}(1:end-4) '_XHz' '.mat'];
            displaytext=['Saving data to:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            save([char(ExportDataDirectory) '\' fname],'DataEventHypnog_Mat','ChannelsList','ColumnHeads');
            displaytext='Finished saving data';
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
        else
            fname=[Filenames{n}(1:end-4) '_XHz' '.mat'];
            displaytext=['User skipped:' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
        end
    catch me
        displaytext=['Error converting: ' fname];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        errorlist = [errorlist n];
    end
end
displaytext='Complete';
disp(displaytext); set(handletext,'String',displaytext); drawnow;
if size(errorlist)>0
    displaytext=['Errors: ' num2str(errorlist)];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = ImportOneToMat(n,Filenames,ChannelNumbers)
system = char(Filenames(n,7));
if strcmp(system,'Alice')
    [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Alice(n,Filenames,ChannelNumbers);
end
if strcmp(system,'ProfusionXML')
    [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = ProfusionXML(n,Filenames,ChannelNumbers);
end
if strcmp(system,'Spike')
    [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Spike(n,Filenames);
end
if strcmp(system,'AliceSleepwareG3')
    [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = AliceSleepwareG3(n,Filenames,ChannelNumbers);
end
if strcmp(system,'GrassTwin')
    [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = GrassTwin(n,Filenames,ChannelNumbers);
end
if strcmp(system,'Remlogic1p1')
    [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Remlogic1p1(n,Filenames,ChannelNumbers);
end
%% Alice
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Alice(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    %EDF Start time
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get recording start time from EDF: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    fid = fopen([directory '\' fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the EDF files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    %F_samp=100;
    
    EDFfilenamedir = [directory '\' fname];
    
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time;
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
    end
    
    
    %Alice Hypnogram
    
    directory = char(Filenames(n,6));
    fname=char(Filenames(n,3));
    displaytext=['Get Hypnogram data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    [num,~,~]=xlsread([directory '\' fname],1,'A:A');
    hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Alice Events
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    Events_fname = [directory '\' fname];
    
    [~,Event_Type,~] = xlsread([Events_fname],'a:a'); % Read the event type column
    Nlines=length(Event_Type);
    Event_Type(1)=[]; % remove the header line of the data
    
    [Event_StartTime,~,~] = xlsread([Events_fname],'c:c'); % Read the column defining the event start time
    Event_StartTime(Event_StartTime<0.5)=Event_StartTime(Event_StartTime<0.5)+1;
    Event_StartTime=Event_StartTime*86400;
    
    [~,~,EventDurationText] = xlsread([Events_fname],['f2:f' num2str(Nlines)]); % Read the column defining the event duration data
    for m=1:length(EventDurationText)
        if length(EventDurationText{m})>1
            Event_Duration{m}=sscanf((EventDurationText{m}),'%f'); %finds first number in text column
        else
            Event_Duration{m}=EventDurationText{m};
        end
    end
    
    %Add Events data in additional colums of the matrix:
    %**************************************************************************
    EventCategories={'Central apnea','Obstructive apnea','Mixed apnea','Hypopnea','Desaturation','µ-arousal'};
    EventCols=[10 11 11 12 13 14];
    
    for m=1:length(Event_Type)
        for i=1:length(EventCategories)
            if strcmp(Event_Type{m},EventCategories{i})==1
                lefti=round((Event_StartTime(m)-StartTime)*F_samp)+1;
                righti=lefti+round(Event_Duration{m}*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
            end
        end
    end
    
    if 0 %check plot
        figure(1)
        plotcols=[2,11,12,3,9,10];
        downsamplefactor=20;
        hold('off')
        count=1;
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
        end
        linkaxes(ax,'x');
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','µ-arousal'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                   3                     2        9           10        11                    12         13       14         4           5       6       7          8      ];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive(/Mixed) 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Xml files from Profusion
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = ProfusionXML(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    %EDF Start time
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get recording start time from EDF: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    fid = fopen([directory '\' fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the EDF files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    
    %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
    
    EDFfilenamedir = [directory '\' fname];
    
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,LabelTemp,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            displaytext=['Found channel labelled:' LabelTemp ' at ' num2str(Fs), ' Hz'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time; 
    
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
        eval(['clear ' Channels{i}]); %save memory
    end
    
    %Xml Hypnogram
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get Hypnogram data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    S=xml2struct([directory '\' fname]);
    
    Nepochs = length(S.CMPStudyConfig.SleepStages.SleepStage);
    clear hypnogram
    hypnogram=NaN*zeros(Nepochs,1);
    for i=1:Nepochs
        hypnogram(i) = str2num(S.CMPStudyConfig.SleepStages.SleepStage{1,i}.Text);
    end %W=0,N1=1,N2=2,N3=3,R=5
    
    if sum(hypnogram==0)==length(hypnogram)
        displaytext=['Warning: Entire hypnogram is wake'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    hypnogram(hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
    hypnogram(hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
    hypnogram=3-hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    clear Pnasal_Time; %save memory
    
    DataEventHypnog_Mat(:,9)=EpochsXHz; 
    clear EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Xml events
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    if ~isfield(S.CMPStudyConfig.ScoredEvents,'ScoredEvent')
        displaytext=['Warning: No scored events'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    Nevents = length(S.CMPStudyConfig.ScoredEvents.ScoredEvent);
 
    clear Events
    for i=1:Nevents
        Events{i}.name = S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Name.Text;
        Events{i}.start = StartTime + str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Start.Text);
        Events{i}.duration = str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Duration.Text);
    end %W=0,N1=1,N2=2,N3=3,R=5
    
    if 0
    for i=1:Nevents
        [num2str(i) ':' Events{i}.name]
    end %
    end
    
    %Add Events data in additional colums of the matrix:
    %**************************************************************************
    EventCategories={'Central Apnea','Obstructive Apnea','Mixed Apnea','Hypopnea','SpO2 desaturation','Arousal','Unsure'}
    
    EventCols=[10 11 15 12 13 14 12]; %Mixed apneas are no longer treated as Obstructive apneas here. Unsure is hypopnea.
    
    %initialize
    DataEventHypnog_Mat(1:size(DataEventHypnog_Mat,1),10:15)=0;
    
    for m=1:Nevents
        for i=1:length(EventCategories)
            strlength=length(EventCategories{i});
            if length(Events{m}.name)>=strlength&&(strcmp(Events{m}.name(1:strlength),EventCategories{i})==1)
                lefti=round((Events{m}.start-StartTime)*F_samp)+1;
                righti=lefti+round((Events{m}.duration)*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
            end
        end
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','µ-arousal','Mixed Apnea'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
    
        
    if 1 %check plot
        figure(1)
        plotcols=[2,3,9,10,11,12,15];
        downsamplefactor=4;
        hold('off')
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
            ylabel(ChannelsList{plotcols(i)});
        end
        linkaxes(ax,'x');
    end
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Spike data
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Spike(n,Filenames)
global handletext F_samp
try
    
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get data from Spike in Matlab format: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    load([directory '\' fname]);
    
    %% Options hardcoded
    MultiplyPositionChannelby=5; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
    
    %% Make channel names uniform: OK
    
    %First name in list is the 'preferred' name.
    %Variables will be renamed to match the preferred name, on a first-come first-served basis
    
    %channel_list = {'Time';'RC';'ABD';'Vflow';'SaO2XHz';'EEG_C3_A2';'PosXHz';'epochsXHz2';'EventsCentralApnea';'EventsObstrApnea';'EventsObstrHyp';'Desats';'EventsAr';'EventsMixedApnea';'EventsCentralHyp';'Pepi';'PCO2';'pO2'};
    clear channelnameoptions
    channelnameoptions.SaO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat'};
    channelnameoptions.Position={'Position','Pos','pos','position'};
    channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal'};
    channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
    channelnameoptions.EEG_O2_A1={'EEG_O2_A1','O2_A1','C4_A2'};
    channelnameoptions.EEG_C3_A2={'EEG','EEG_C3_A2','C3_A2'};
    channelnameoptions.RC={'Thorax','RC','Chest','CHEST','Belt2'};
    channelnameoptions.ABD={'Abdomen','ABD','Abdom','ABDM','Belt1'};
    channelnameoptions.EKG={'EKG','ECG'};
    channelnameoptions.Vflow={'Pnasal','PNasal','Pmask','PMask','Vflow'}; %'Vflow'
    channelnamestemp=fieldnames(channelnameoptions);
    
    for i=1:length(channelnamestemp)
        temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
        for n=1:length(temp)
            if n==1&&exist(char(temp(n)),'var'); %name is already first in the list - stop
                break
            elseif exist(char(temp(n)),'var'); %rename to that which is the first in the list
                eval([char(temp(1)) '=' char(temp(n))]);
                eval(['clear ' char(temp(n))]);
                break
            end
        end
    end
    
    %Delete these channels if they exist
    unusedchannellist={'ECG','EEG_O2_A1','Therm','PCO2','pO2','Snore','snore','F3_A2','F4_A1','O1_A2','Pepi','EPI','Epi','ROC','EOG_R','LOC','EOG_L','EMG','EMGchin','EMG_DIA','EMG_GG','RLeg','LLeg','Keyboard','blank','mergeinfo'}; 
    for i=1:length(unusedchannellist)
        eval(['clear ' unusedchannellist{i} ';']);
    end
    
    %channel_list = {'Time';'RC';'ABD';'Vflow';'SaO2';'EEG_C3_A2';'Position';'Epochs';'EventsCentralApnea';'EventsObstrApnea';'EventsObstrHyp';'Desats';'EventsAr';'EventsMixedApnea';'EventsCentralHyp';'Pepi';'PCO2';'pO2';'ECG';'EMGchin'};
    
    Channels={'Pnasal','Thorax','Abdomen','SaO2','EEG','Position','CPAP'};
    %F_samp=100; 
    dt=1/F_samp;
    
    for i=1:length(Channels)
        displaytext=['Collecting channel: ' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval(['Fs=(1/' Channels{i} '.interval);']);
            if round(Fs)~=F_samp
                displaytext=['Resampling: ' Channels{i} ' from ' num2str(round(Fs)) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} '.values = resample(' Channels{i} '.values,F_samp,round(Fs));']); %only works if Fsamp / Fs are integer multiples of each other
                eval([Channels{i} '.interval = 1/F_samp;']); %only works if Fsamp / Fs are integer multiples of each other
                eval([Channels{i} '.length = length(' Channels{i} '.values);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel: ' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values = -1+ 0*Pnasal.values;']); %assumes we at least have a Pnasal channel
            eval([Channels{i} '.interval = 1/F_samp;']); 
            eval([Channels{i} '.length = length(' Channels{i} '.values);']); 
        end
    end
    
    %% Get timing, epoch and event data
    
    %get time data
    timeXHz=(Pnasal.start:(1/F_samp):(Pnasal.start+(length(Pnasal.values)/F_samp)-1/F_samp))';
    N=length(timeXHz);
    
    %reserve memory for the large matrix
    DataEventHypnog_Mat=zeros(N,16);
        DataEventHypnog_Mat(:,1)=timeXHz;
        clear timeXHz;
        
%% New: fix channel lengths if not = N (not the same as Pnasal) -- incorporate this into Alice and others
    for i=1:length(Channels)
        eval(['ChannelN(i)=' Channels{i} '.length;']);
        if ChannelN(i)>N
            displaytext=['Length of ' Channels{i} ' channel is being altered from ' num2str(ChannelN(i)) ' to ' num2str(N)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values(N+1:end)=[];']);
            eval([Channels{i} '.length=N;']);
        end
        if ChannelN(i)<N
            displaytext=['Length of ' Channels{i} ' channel is being altered from ' num2str(ChannelN(i)) ' to ' num2str(N)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values((ChannelN(i)+1):N)=' Channels{i} '.values(end);']);
            eval([Channels{i} '.length=N;']);
        end
    end
    
    %% Epoch info
    if exist('Epochs','var')
        displaytext=['Get Hypnogram data:' fname];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
        tempdiff=round(double(Epochs.times));
        tempdiff2=diff(tempdiff);
        CountEpochsNot30s=sum(tempdiff2(1:end-1)~=30);
        if CountEpochsNot30s==0
            epochsXHz.values=interp1((Epochs.times+15)',double(Epochs.codes(:,1)'),DataEventHypnog_Mat(:,1),'nearest','extrap');
        else
            epochsXHz.values=zeros(1,N);
            for i=1:length(Epochs.times)
                if i<length(Epochs.times)
                    epochsXHz.values(DataEventHypnog_Mat(:,1)>=Epochs.times(i)&DataEventHypnog_Mat(:,1)<Epochs.times(i+1))=double(Epochs.codes(i,1));%(Epochs.codes(i,1));
                elseif i==length(Epochs.times)
                    %epochsXHz.values(timeXHz>=Epochs.times(i))=double(Epochs.codes(i));
                end
            end
        end
        
        recode=1
        if recode
            %Spike code: W=0,N1=1,N2=2,N3=3,R=5,?=8
            %recode Epochs to match Phil's naming: N3=0,N2=1,N1=2,R=3,W=4
            Epochs.values=3-epochsXHz.values;
            Epochs.values(epochsXHz.values==0)=4;
            Epochs.values(epochsXHz.values==5)=3;
            Epochs.values(epochsXHz.values==8)=NaN;
        end
        clear epochsXHz;
    end
    
    
    %% Events info
    
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    if exist('Evts','var'), New_Evts=Evts; end
    if exist('Evts2','var'), NewEvts2=Evts2; end
    if exist('New_Evts2','var'), NewEvts2=New_Evts2; end
    
    if exist('NewEvts2','var')&&exist('New_Evts','var')
        %combine events lists
        New_Evts_backup1=New_Evts;
        New_Evts.times=[New_Evts.times;NewEvts2.times];
        New_Evts.codes=[New_Evts.codes;NewEvts2.codes];
        New_Evts.text=[New_Evts.text;NewEvts2.text];
        New_Evts.length=New_Evts.length+NewEvts2.length;
    end
    
    if New_Evts.codes(1,1)==0
        New_Evts.codes(1,:)=[];
        New_Evts.times(1)=[];
        New_Evts.text(1,:)=[];
        New_Evts.length=New_Evts.length-1;
        %'removed_one'
    end
    
    ii=1; %loop cycles multiple times, removing any double rows of zeros.
    while ii<length(New_Evts.codes)
        ii=1;
        for i=2:length(New_Evts.times)
            evennumber=1-mod(i,2);
            oddnumber=1-evennumber;
            if 1
                if oddnumber&&New_Evts.codes(i)==0&&New_Evts.codes(i-1)==0
                    New_Evts.codes(i,:)=[];
                    New_Evts.times(i)=[];
                    New_Evts.text(i,:)=[];
                    New_Evts.length=New_Evts.length-1;
                    %'removed_one double row of zeros'
                    break
                end
            end
            %     if evennumber&New_Evts.codes(i)>0
            %         'saysomething'
            %         break
            %     end
            ii=ii+1;
        end
    end
    
    % insert rows of missing zeros.
    for i=2:length(New_Evts.times)
        evennumber=1-mod(i,2);
        oddnumber=1-evennumber;
        if evennumber&&New_Evts.codes(i)>0&&New_Evts.codes(i-1)>0
            if 1
                New_Evts.codes=[New_Evts.codes(1:(i-1),:);[0 0 0 0];New_Evts.codes(i:end,:)];
                New_Evts.times=[New_Evts.times(1:i);New_Evts.times(i:end)];
                New_Evts.text=[New_Evts.text(1:(i-1),:);['                                            '];New_Evts.text(i:end,:)];
                New_Evts.length=New_Evts.length+1;
                %'added_row_of_zeros'
            end
        end
    end
    
    %% Make list of event types
    
    clear EventTypeList_text EventTypeList_code tempi
    EventTypeList_text=[];
    for i=1:2:size(New_Evts.codes,1)
        codematch=0;
        for j=1:size(EventTypeList_text,1)
            if strcmp(char(EventTypeList_text(j,:)),New_Evts.text(i,:))
                %'found code match'
                codematch=1;
                break
            end
        end
        if codematch==0
            tempi=size(EventTypeList_text,1);
            EventTypeList_code(tempi+1)=New_Evts.codes(i);
            EventTypeList_text(tempi+1,:)=New_Evts.text(i,:);
        end
    end
    EventTypeList_code=EventTypeList_code'
    char(EventTypeList_text)
    
    %Make list fit the following categories:
    %1 AR
    %2 ApO | Ap-O
    %3 ApC | Ap-C
    %4 HypO | H (Y)
    %5 M (Y)
    %6 HypC
    %7 HypO2% (Y)
    %8 HypC2%
    %9 HypOx (Y)
    %10 HypCx
    
    New_Evts_backup2=New_Evts;
    
    for i=1:size(New_Evts.codes,1)
        if strncmp(New_Evts.text(i,:),'AR',2)
            New_Evts.codes(i,1)=1;
        elseif strncmp(New_Evts.text(i,:),'ApO',3)
            New_Evts.codes(i,1)=2;
        elseif strncmp(New_Evts.text(i,:),'Ap-O',4)
            New_Evts.codes(i,1)=2;
        elseif strncmp(New_Evts.text(i,:),'ap-O',4)
            New_Evts.codes(i,1)=2;
        elseif strncmp(New_Evts.text(i,:),'ApC',3)
            New_Evts.codes(i,1)=3;
        elseif strncmp(New_Evts.text(i,:),'Ap-C',4)
            New_Evts.codes(i,1)=3;
        elseif strncmp(New_Evts.text(i,:),'M',1)
            New_Evts.codes(i,1)=5;
        elseif strncmp(New_Evts.text(i,:),'HypO2',5)
            New_Evts.codes(i,1)=7;
        elseif strncmp(New_Evts.text(i,:),'HypC2',5)
            New_Evts.codes(i,1)=8;
        elseif strncmp(New_Evts.text(i,:),'HypOx',5)
            New_Evts.codes(i,1)=9;
        elseif strncmp(New_Evts.text(i,:),'HypCx',5)
            New_Evts.codes(i,1)=10;
        elseif strncmp(New_Evts.text(i,:),'HypO',4)
            New_Evts.codes(i,1)=4;
        elseif strncmp(New_Evts.text(i,:),'HypC',4)
            New_Evts.codes(i,1)=6;
        elseif strncmp(New_Evts.text(i,:),'Hyp-C',5)
            New_Evts.codes(i,1)=6;
        elseif strncmp(New_Evts.text(i,:),'H',1)
            New_Evts.codes(i,1)=4;
        end
    end
    
    %%%Add this in%%%
    % for x=1:2:length(New_Evts.times)
    %     New_Evts.duration((x-1)/2+1)=New_Evts.times(x+1)-New_Evts.times(x);
    % end
    
    %% Make events in continuous time
    % if exist('New_Evts','var')
    
    %channel_list = {... 'EventsCentralApnea';'EventsObstrApnea';'EventsObstrHyp';'Desats';'EventsAr';'EventsMixedApnea';'EventsCentralHyp'};
    EventsAr.values=zeros(N,1); EventsAr.interval=dt;
    EventsObstrApnea.values=zeros(N,1); EventsObstrApnea.interval=dt;
    EventsCentralApnea.values=zeros(N,1); EventsCentralApnea.interval=dt;
    EventsObstrHyp.values=zeros(N,1); EventsObstrHyp.interval=dt;
    EventsMixedApnea.values=zeros(N,1); EventsMixedApnea.interval=dt;
    EventsCentralHyp.values=zeros(N,1); EventsCentralHyp.interval=dt;
    
    for x=1:2:length(New_Evts.times)
        lefti=round((New_Evts.times(x)-DataEventHypnog_Mat(1,1))*F_samp)+1;
        righti=lefti+round((New_Evts.times(x+1)-New_Evts.times(x))*F_samp);
        if New_Evts.codes(x)==1
            EventsAr.values(lefti:righti)=New_Evts.codes(x);
        end
        if New_Evts.codes(x)==2
            EventsObstrApnea.values(lefti:righti)=1;
        end
        if New_Evts.codes(x)==3
            EventsCentralApnea.values(lefti:righti)=1;
        end
        if New_Evts.codes(x)==4
            EventsObstrHyp.values(lefti:righti)=1;
        end
        if New_Evts.codes(x)==5
            EventsMixedApnea.values(lefti:righti)=1;
        end
        if New_Evts.codes(x)==6
            EventsCentralHyp.values(lefti:righti)=1;
        end
    end
     
    % Plot events
    if 0
        figure(1);
        plot(DataEventHypnog_Mat(:,1),EventsAr.values)%,New_Evts.times,New_Evts.codes(:,1),'.'); ylabel('EventsIndex'); end
        set(gca,'FontName','Arial Narrow','FontSize',10,'box','off','YTickLabel',{'AR','OA','CA','OH','M','CH','OH2','CHx','OHx','CHx'},'YTick',[1 2 3 4 5 6 7 8 9 10],'XTick',[],'Xcolor',[1 1 1]);
        ylabel('Events/Arousal');
    end
    
    %% Position data at XHz.
    
    clear Mexception
    try %if channel doesn't exist: program will continue
        channel_name = 'Position';
        Position.values=Position.values*MultiplyPositionChannelby;       
        displaytext='Scaled position data';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;    
    catch Mexception
        displaytext='No position data, looking for file';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
    
    % Fix Pos data if necessary
    textfilename=[directory '\' fname(1:(length(fname)-4)) '_pos.txt'];
    if exist(textfilename,'file')==2
        displaytext='Found position file';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        
        %PosXHz_original=Position.values;
        [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
        for i=1:(length(col1)-1)
            Position.values(DataEventHypnog_Mat(:,1)>=col1(i)&DataEventHypnog_Mat(:,1)<col1(i+1))=col2(i);
        end
        Position.values(DataEventHypnog_Mat(:,1)>=col1(end))=col2(end);
    end
    %channel_exist=sum(abs(Position.values))>0;
    
    if 0
        figure(1), plot(DataEventHypnog_Mat(:,1),Position.values);
        set(gca,'FontName','Arial Narrow','FontSize',10,'box','off','YTickLabel',{'Up','Left','Supine','Prone','Right'},'YTick',[0 1 2 3 4],'XTick',[],'Xcolor',[1 1 1],'Ylim',[0 5]);
        ylabel('Position');
    end
    %clear PosXHz_original

    
    %% Channel data at XHz.
        %DataEventHypnog_Mat=[Flow_Aux_Time.' RIP_Thorax RIP_Abdo Flow_Aux SpO2 C3A2 Position CPAP];
        ColumnHeads=[1      2            4      8        9               10            11         12       13         3           5       6           7          14             15           16            17       18    19];
        %            [1=Time 2=RIP_Thorax 3=Vflow 4=Hypnog 5=CentralApnea 6=Obstr.Apnea 7=ObsHypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=EEGC3A2 13=Position 14=CPAP 15=Mixed apnea 16=CentrHypop 17=Pepi/Pes 18=PCO2 19=PO2]
        
        
        %ColumnHeads=[1      2            4      8        9         10            11         12       13         3           5       6       7];
        %            [1=Time 2=RIP_Thorax 3=Flow 4=Hypnog 5=Central 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]
        
        ChannelsList = {'Time';'Thorax';'Abdomen';'Pnasal';'SaO2';'EEG';'Position';'Epochs';'EventsCentralApnea';'EventsObstrApnea';'EventsObstrHyp';'Desats';'EventsAr';'CPAP';'EventsMixedApnea';'EventsCentralHyp'};
        
        displaytext='Combining data into large matrix';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        for x=2:length(ChannelsList)
            try %if channel doesn't exist: program will continue but an error message will show
                eval(['DataEventHypnog_Mat(:,x)=' ChannelsList{x} '.values;']);               
            catch me
                displaytext=['No channel found: ' me.message];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
            end
        end
        
        if 0
            plotchannels=[2 3 4 6]
            figure(1);
            for i=1:length(plotchannels)
            ax(i)=subplot(length(plotchannels),1,i); 
            plot(DataEventHypnog_Mat(:,1),DataEventHypnog_Mat(:,plotchannels(i)))
            end
            linkaxes(ax,'x');
        end
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Alice Sleepware G3
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = AliceSleepwareG3(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get rml data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    S=xml2struct([directory '\' fname]);
    
    displaytext=['Get recording start time from rml: ' directory];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    StartText = S.PatientStudy.Acquisition.Sessions.Session.RecordingStart.Text(12:19);
    StartTime = mod(datenum(StartText,'HH:MM:SS'),1)*86400
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the .csv files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    NChannels = length(Channels);
    %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
    %EDF Start time -- have to import codes
    directory = char(Filenames(n,4));
    
    fname=char(Filenames(n,1));
    EDFfilenamedir = [directory '\' fname];
    
    %M = csvread(filename) 
    if 1
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    end
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time; 
    %clear Pnasal_Time; %save memory
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
        eval(['clear ' Channels{i}]); %save memory
    end
    
    %Xml Hypnogram
    displaytext=['Get Hypnogram data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
   
    StageTextOptions = {'Wake','NonREM1','NonREM2','NonREM3','REM','NotScored'};
    StageTextCodes = [4,2,1,0,3,-1];
    Ntransitions = length(S.PatientStudy.ScoringData.StagingData.UserStaging.NeuroAdultAASMStaging.Stage);
    clear Transitions
    for i=1:Ntransitions
        Transitions{i}.Text = S.PatientStudy.ScoringData.StagingData.UserStaging.NeuroAdultAASMStaging.Stage{i}.Attributes.Type;
        Transitions{i}.Time = str2num(S.PatientStudy.ScoringData.StagingData.UserStaging.NeuroAdultAASMStaging.Stage{i}.Attributes.Start)/30 +1;
        Transitions{i}.Code = -1; %Default is unknown
    end 
    for i=1:Ntransitions
        for m=1:length(StageTextOptions)
            if strcmp(StageTextOptions{m},Transitions{i}.Text)==1
                %[StageTextOptions{m} ' ' num2str(StageTextCodes(m))]
                Transitions{i}.Code = StageTextCodes(m);
            end
        end
    end
    
    RecordingDuration = str2num(S.PatientStudy.Acquisition.Sessions.Session.Duration.Text);
    Nepochs = ceil(RecordingDuration/30);
    
    clear hypnogram
    hypnogram=NaN*zeros(Nepochs,1);
    for i=1:Ntransitions
        indexi = Transitions{i}.Time;
        hypnogram(indexi:end) = Transitions{i}.Code;
    end 
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz; 
    clear EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Xml/rml events [tested and working]
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    
    TempEvents = S.PatientStudy.ScoringData.Events.Event;
    EventTypesOfInterest = {'ObstructiveApnea','Hypopnea','CentralApnea','MixedApnea','Arousal','RelativeDesaturation'};

    Nevents = length(S.PatientStudy.ScoringData.Events.Event);
    clear Events
    for i=1:Nevents
        Events{i}.name = S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Type;
        Events{i}.start = StartTime + str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Start);
        Events{i}.duration = str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Duration);
    end %W=0,N1=1,N2=2,N3=3,R=5
    
    %Add Events data in additional colums of the matrix:
    %**************************************************************************
    EventCategories={'CentralApnea','ObstructiveApnea','MixedApnea','Hypopnea','RelativeDesaturation','Arousal'}
    
    EventCols=[10 11 15 12 13 14]; %Mixed apneas are no longer treated as Obstructive apneas here.
    
    for m=1:Nevents
        for i=1:length(EventCategories)
            if strcmp(Events{m}.name,EventCategories{i})==1
                %Events{m}.name
                lefti=round((Events{m}.start-StartTime)*F_samp)+1;
                righti=lefti+round((Events{m}.duration)*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
            end
        end
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','µ-arousal'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
    
        
    if 1 %check plot
        figure(1)
        plotcols=[2,11,12,3,9,7];
        downsamplefactor=4;
        hold('off')
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
            ylabel(ChannelsList{plotcols(i)});
        end
        linkaxes(ax,'x');
    end
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% GrassTwin [needs checking]
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = GrassTwin(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get xls data: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    [~,xlstxt,~] = xlsread([directory '\' fname],1,'C:D')
    
   
    displaytext=['Get recording start time from xlstxt'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    StartText = xlstxt{1,1};
    StartTime = mod(datenum(StartText,'HH:MM:SS'),1)*86400
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the .csv files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    NChannels = length(Channels);
    %F_samp=100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normally 100
    %EDF Start time -- have to import codes
    directory = char(Filenames(n,4));
    
    fname=char(Filenames(n,1));
    EDFfilenamedir = [directory '\' fname];
    
    %M = csvread(filename) 
    if 1
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    end
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time; 
    %clear Pnasal_Time; %save memory
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
        eval(['clear ' Channels{i}]); %save memory
    end
    
    %Xls Hypnogram
    displaytext=['Get Hypnogram data from xlstxt'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
   
    StageTextOptions = {'Stage - W','Stage - R','Stage - N1','Stage - N2','Stage - N3','Stage - No Stage'};
    StageTextCodes = [4,3,2,1,0,-1]; %Terrill code: W=4,R=3,N1=2,N2=1,N3=0 [unknown=-1]

    xlstxthyponly = xlstxt;
    temp = 'Stage - ';
    for i=size(xlstxthyponly,1):-1:1
        temphyp=xlstxthyponly{i,2};
        if length(temphyp)<8||~strcmp(temphyp(1:8),'Stage - ')
            xlstxthyponly(i,:)=[];
        end
    end
    Ntransitions = size(xlstxthyponly,1);
    LastEpochStartTime = mod(datenum(xlstxthyponly{Ntransitions,1},'HH:MM:SS'),1)*86400
    if LastEpochStartTime<43200, LastEpochStartTime=LastEpochStartTime+86400; end
    Nepochs = round((LastEpochStartTime-StartTime)/30+1);
    
    clear HypTimes
    for i=1:Ntransitions
    HypTimes(i,:) = mod(datenum(xlstxthyponly{i,1},'HH:MM:SS'),1)*86400
    if HypTimes(i)<43200, HypTimes(i)=HypTimes(i)+86400; end
    end
    
    clear HypCodes
    for i=1:Ntransitions
        for m=1:length(StageTextOptions)
            if strcmp(xlstxthyponly{i,2},StageTextOptions{m})
            HypCodes(i,:) = StageTextCodes(m);
            end
        end
    end
    
    HypEpochIndex = round((HypTimes-StartTime)/30+1);
    %RecordingDuration = ... **might be more robust to use the EDF data for this in future
    %Nepochs = ceil(RecordingDuration/30);
    
    clear hypnogram
    hypnogram=NaN*zeros(Nepochs,1);
    for i=1:Ntransitions %write from current epoch to end file with new stage code
        indexi = round(HypEpochIndex(i));
        hypnogram(indexi:end) = HypCodes(i);
    end 
    
    hypnogram_t = StartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz; 
    clear EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Xls events 
    displaytext=['Get Events data from xlstxt'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    
    %xlstxt
    EventTypesOfInterest = {'Central Apnea','Obstructive Apnea','Mixed Apnea','Hypopnea','RERA','Desaturation','Arousal'};
    EventCodes = [10 11 15 12 12 13 14]; %RERAs are treated as hypopneas here
    %Nevents = length(S.PatientStudy.ScoringData.Events.Event);
    index = 1;
    clear Events
    for i=1:size(xlstxt,1)
        for m=1:length(EventTypesOfInterest)
            if ~isempty(findstr(xlstxt{i,2},EventTypesOfInterest{m}))
                Events.name{index}=EventTypesOfInterest{m};
                Events.code(index)=EventCodes(m);
                Events.row(index)=i;
                timetemp = mod(datenum(xlstxt{i,1},'HH:MM:SS'),1)*86400;
                    if timetemp<43200, timetemp=timetemp+86400; end
                    Events.start(index)=timetemp;
                tempk1 = findstr(xlstxt{i,2},'Dur:') + 5;
                    tempk2 = findstr(xlstxt{i,2},'sec') - 2;
                    Events.duration(index) = str2num(xlstxt{i,2}(tempk1:tempk2));
                index=index+1;    
            end
            continue
        end
    end 
    Nevents = length(Events.code);
    
    %Add Events data in additional colums of the matrix:   
    for i=1:Nevents
          lefti=round((Events.start(i)-StartTime)*F_samp)+1;
          righti=lefti+round((Events.duration(i))*F_samp);
          DataEventHypnog_Mat(lefti:righti,Events.code(i))=1;
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','µ-arousal'}];
    
    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
    
        
    if 1 %check plot
        figure(1)
        plotcols=[2,11,12,3,9,7];
        downsamplefactor=4;
        hold('off')
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
            ylabel(ChannelsList{plotcols(i)});
        end
        linkaxes(ax,'x');
    end
    
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Alice
function [DataEventHypnog_Mat,ChannelsList,ColumnHeads] = Remlogic1p1(n,Filenames,ChannelNumbers)
global handletext F_samp

try
    %EDF Start time
    directory = char(Filenames(n,4));
    fname=char(Filenames(n,1));
    
    displaytext=['Get recording start time from EDF: ' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    fid = fopen([directory '\' fname],'r');
    fseek(fid,176,-1);
    StartTimeText = char(fread(fid,8,'char')');
    fclose(fid); % Close file
    StartTime = mod(datenum(StartTimeText,'HH.MM.SS'),1)*86400;
    if StartTime<43200, StartTime=StartTime+86400; end
    
    % This section imports respiratory flow and other polysomnography signal data from the EDF files
    Channels={'Pnasal','Thorax','Abdomen','SpO2','EEG','Position','CPAP'};
    %F_samp=100;
    
    EDFfilenamedir = [directory '\' fname];
    
    for i=1:length(Channels)
        displaytext=['Collecting channel:' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            eval(['[' Channels{i} ',Fs,~,~,~,~,~,~,~] = readedfrev3(EDFfilenamedir,ChannelNumbers(n,i)-1,0,Inf);']) %the '-1' is because EDF channel numbers start from 0.
            %e.g. [Pnasal,Fs,~,~,~,~,~,~,~] = readedfrev3([directory '\' fname],ChannelNumbers(i)-1,0,Inf); %Extract 100Hz Flow data
            if Fs~=F_samp
                displaytext=['Resampling:' Channels{i} ' from ' num2str(Fs) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} ' = resample(' Channels{i} ',F_samp,Fs);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel:' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} ' = -1+ 0*Pnasal;']); %assumes we at least have a Pnasal channel
        end
    end
    
    N_flow=length(Pnasal);
    Pnasal_Time=(StartTime:(1/F_samp):StartTime+(N_flow-1)*(1/F_samp))'; % This is the time vector associated with the 100Hz Flow data.
    
    for i=1:length(Channels)
        N_Channels(i)=length(eval(Channels{i}));
    end
    
    
    %Combine data into a single matrix
    displaytext='Combining data into a single matrix';
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    DataEventHypnog_Mat=Pnasal_Time;
    for i=1:length(Channels)
        eval(['DataEventHypnog_Mat=[DataEventHypnog_Mat ' Channels{i} '];']);%
    end
       
    %Hypnogram
    directory = char(Filenames(n,6));
    fname=char(Filenames(n,3));
    displaytext=['Get Hypnogram data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
    clear fileinfo
        fid = fopen([directory '\' fname]);
        i=1;
        while 1
            fileinfo{i,:} = fgetl(fid);  %read in the data
            if fileinfo{i,:}==-1
                fileinfo(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        starttext = 'Sleep Stage';
         %%%up to here
        for i=1:length(fileinfo)
           if ~isempty(findstr(fileinfo{i,:},starttext)) 
               break
           end
        end
        startline = i+1;
        for i=1:length(fileinfo)
           if ~isempty(findstr(fileinfo{i,:},starttext)) 
               break
           end
        end
        sleeptable = fileinfo(startline:end);
       
        clear SleepCodesTxt
    for i=1:length(sleeptable);
        SleepCodesTxt{i} = sleeptable{i}(1:2);
    end
        
    StageTextOptions = {'W','R','N1','N2','N3'};
    StageTextCodes = [4,3,2,1,0];
    
    clear hypnogram
    for i=1:length(SleepCodesTxt)
        for m=1:length(StageTextOptions)
            if findstr(SleepCodesTxt{i},StageTextOptions{m})
                   hypnogram(i)=StageTextCodes(m);
                   break
            end
        end
    end
    
    HypStartTimeTxt = sleeptable{1}(11:18);
    
    timestrrangetemp=findstr(sleeptable{1},':');
        timestrrange = timestrrangetemp(1)-2:timestrrangetemp(2)+2;
        HypStartTimeTxt = sleeptable{1}(timestrrange);
        
    HypStartTime = mod(datenum(HypStartTimeTxt,'HH:MM:SS'),1)*86400;
    if HypStartTime<43200, HypStartTime=HypStartTime+86400; end
    
    hypnogram_t = HypStartTime+(0:30:(30*(length(hypnogram)-1)))';
    EpochsXHz = interp1(hypnogram_t+15,hypnogram,Pnasal_Time,'nearest','extrap');
    
    DataEventHypnog_Mat(:,9)=EpochsXHz;
    
    %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
    
    %Events
    directory = char(Filenames(n,5));
    fname = char(Filenames(n,2));
    displaytext=['Get Events data:' fname];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    Events_fname = [directory '\' fname];
    
    clear fileinfo
        fid = fopen(Events_fname);
        i=1;
        while 1
            fileinfo{i,:} = fgetl(fid);  %read in the data
            if fileinfo{i,:}==-1
                fileinfo(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        starttext = 'Sleep Stage';
        for i=1:length(fileinfo)
           if ~isempty(findstr(fileinfo{i,:},starttext)) 
               break
           end
        end
        startline = i+1;
        for i=1:length(fileinfo)
           if ~isempty(findstr(fileinfo{i,:},starttext)) 
               break
           end
        end
        EventsTable = fileinfo(startline:end);
       

    %up to here
    clear Event_StartTime Event_StartTimeTxt
    for i=1:length(EventsTable)
        timestrrangetemp=findstr(EventsTable{i},':');
        timestrrange = timestrrangetemp(1)-2:timestrrangetemp(2)+2;
        Event_StartTimeTxt{i} = EventsTable{i}(timestrrange);
        Event_StartTime(i,:) = mod(datenum(Event_StartTimeTxt{i},'HH:MM:SS'),1)*86400;
            if Event_StartTime(i)<43200, Event_StartTime(i)=Event_StartTime(i)+86400; end
    end
    
    clear Event_DurationText Event_Duration
    for i=1:length(EventsTable)
        Event_DurationText{i}=EventsTable{i}(end-1:end);
        Event_Duration(i,:)=str2num(Event_DurationText{i});
    end
    
    clear Event_StartTime Event_StartTimeTxt
    for i=1:length(EventsTable)
        timestrrangetemp=findstr(EventsTable{i},':');
        timestrrange = timestrrangetemp(1)-2:timestrrangetemp(2)+2;
        Event_StartTimeTxt{i} = EventsTable{i}(timestrrange);
        Event_StartTime(i,:) = mod(datenum(Event_StartTimeTxt{i},'HH:MM:SS'),1)*86400;
            if Event_StartTime(i)<43200, Event_StartTime(i)=Event_StartTime(i)+86400; end
    end
    %up to here
        EventCategories={'APNEA-CENTRAL','APNEA-OBSTRUCTIVE','APNEA-MIXED','HYPOPNEA','AROUSAL'}; %note there is no desat column but an error will occur if it is missing
        EventCols=[10 11 15 12 14];
       
    clear Event_Type
    for m=1:length(EventsTable)
        for i=length(EventCategories):-1:1
            temp = findstr(EventsTable{m},EventCategories{i});
            if ~isempty(temp)
                Event_Type{m,:}=EventsTable{m}(temp:temp+length(EventCategories{i})-1);
                if EventCols(i)==14&&Event_Duration(m)<3
                    Event_Duration(m)==3; %force arousal duration to be at least 3 s
                end
                lefti=round((Event_StartTime(m)-StartTime)*F_samp)+1;
                righti=lefti+round(Event_Duration(m)*F_samp);
                DataEventHypnog_Mat(lefti:righti,EventCols(i))=1;
                break
            end
        end
    end
    
    if 0 %check plot
        figure(1)
        plotcols=[2,11,12,14,3,9];
        downsamplefactor=5;
        hold('off')
        count=1;
        clear ax
        for i=1:length(plotcols)
            ax(i)=subplot(length(plotcols),1,i);
            plot(downsample(DataEventHypnog_Mat(:,1),downsamplefactor),downsample(DataEventHypnog_Mat(:,plotcols(i)),downsamplefactor));
        end
        linkaxes(ax,'x');
    end
    
    ChannelsList = [{'Time'},Channels,{'Hypnogram','Central apnea','Obstructive apnea','Hypopnea','Desaturation','µ-arousal'}];

    %Set Columnheads, a numerical list of channels, which reorders ChannelsList in the order expected by Analysis.m
    ColumnHeads=[1                       3      2        9           10               11            12         13       14         4           5       6       7           8       15];
    %[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position 14=CPAP 15=Mixed apnea ]
   
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end
