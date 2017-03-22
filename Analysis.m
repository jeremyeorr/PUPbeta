function Analysis()

global AnalyzeDataSpreadsheet handletext

displaytext='Starting up analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;

if isempty(AnalyzeDataSpreadsheet)
    AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
end
[~,txt,raw] = xlsread(AnalyzeDataSpreadsheet,2,'B3:C22');

global settings


%globals include: savename Pnasaldownisinsp LGfromFlowVersion saveplots sqrt_scaling plotfigure usescoredcentralapneas eventsarebreathsfullywithinmargins havescoredcentralhypops manualscoringtouchups maxdelaybreaths Fs exportresultstoxls

settings.savename = char(txt(1,2));
settings.OutputDataDirectory = char(txt(18,2));
try
    load([settings.OutputDataDirectory '\' settings.savename],'settings','AnalysisIndex','LGplusinfo','EventsInfo','LG_QualityInfo','DataOut','ArousalDat','fitQual','SleepData');
    displaytext=['Loading saved data from ' settings.savename];
catch me
    displaytext=['No saved data file called' settings.savename];
end
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    
settings.Pnasaldownisinsp=cell2mat(raw(2,2));
settings.LGfromFlowVersion = char(txt(3,2));
settings.sqrt_scaling=cell2mat(raw(4,2));
settings.saveplots=cell2mat(raw(5,2));
settings.plotfigure=cell2mat(raw(6,2));
settings.usescoredcentralapneas=cell2mat(raw(7,2));
settings.eventsarebreathsfullywithinmargins=cell2mat(raw(8,2));
settings.havescoredcentralhypops=cell2mat(raw(9,2));
settings.manualscoringtouchups=cell2mat(raw(10,2));
settings.maxdelaybreaths=cell2mat(raw(11,2));
settings.exportresultstoxls=cell2mat(raw(15,2));
settings.handlemixedeventsseparately=cell2mat(raw(16,2));
settings.longestwakeduration=cell2mat(raw(17,2));

settings.Fs=cell2mat(raw(13,2));
settings.ignoreCPAPdata=cell2mat(raw(14,2));
settings.locttot=0;

[num,patients] = xlsread(AnalyzeDataSpreadsheet,1,'B3:E5003');
analyzelist = num(:,2);
invertflowlist = num(:,1);
settings.invertflowlist=find(invertflowlist'==1);
settings.windowlength=cell2mat(raw(12,2));
settings.Fs=cell2mat(raw(13,2));

settings.ARmodel=cell2mat(raw(19,2));
settings.AnalyzeNREMonly=cell2mat(raw(20,2));

%%
global n
settings.rerunspecificwindows=[]; %leave empty, work is needed so that if this is saved, other saved data are preserved

WindowDuration=settings.windowlength*60; % Duration of analysis window
settings.WindowStep=120;
for n=1:size(patients,1); %%%%%%%%%%%%%%%%%%%
    if analyzelist(n)==0
        displaytext=['Skipping: n=' num2str(n) ', ' char(patients(n,1))];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        continue
    end
    displaytext=['Patient ' num2str(n) ': ' char(patients(n,1))];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    %Events_fname=O2_files{n,(m-1)*9+2}
    clear ColumnHeads
    
    try
        directoryn=char(patients(n,2));
        if directoryn(end)~='\'
            directoryn=[directoryn '\'];
        end
        MATfilename=[directoryn char(patients(n,1))];
        load(MATfilename);
        
        if ~exist('ColumnHeads','var')
            ColumnHeads=[1    2          4    9      10      11          12       13     14       3        5    6    7          8      ];
            %[1=Time 2=RIP_Thorax 3=Flow 4=Hypnog 5=Central 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]
        end
        
        
        %clear AnalysisIndex LGplusinfo EventsInfo LG_QualityInfo DataOut ArousalDat fitQual
        [AnalysisIndex{n},LGplusinfo{n},EventsInfo{n},LG_QualityInfo{n},DataOut{n},ArousalDat{n},fitQual{n},SleepData{n}]...
            = WindowSelectAndRun(DataEventHypnog_Mat,settings.WindowStep,settings.ARmodel,WindowDuration,ColumnHeads) %pi/16
        
        if settings.exportresultstoxls
            %% Write to xls
            copyfile('ResultsTemplate.xlsx',[settings.OutputDataDirectory '\' char(patients(n,1)) '_results.xlsx'])
            
            displaytext=['Writing to ' char(patients(n,1)) '_results.xlsx'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            
            cell0='H12'; %for filename
            temp_fitQual=[];
            for i=1:length(fitQual{n})
                if ~isempty(fitQual{n}{i})
                    temp_fitQual(i,1:3)=fitQual{n}{i};
                else
                    temp_fitQual(i,1:3)=[NaN NaN NaN];
                end
            end
            try
            exportdatamatrix=[(1:(size(LGplusinfo{n},1)))' LGplusinfo{n}(:,[1 7 8 2 3 9 10 12 13 6 5]) temp_fitQual(:,2) LG_QualityInfo{n}(:,[1 2 3 5 8 11 12]) SleepData{n}]; %will cause an error if you don't run every window...
            xlswrite([settings.OutputDataDirectory '\' char(patients(n,1)) '_results.xlsx'],exportdatamatrix,1,cell0);
            catch me
                disp(me.message)
            end
        end
        
        displaytext=['Save analysis data to ' settings.savename '.mat'];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        save([settings.OutputDataDirectory '\' settings.savename],'settings','AnalysisIndex','LGplusinfo','EventsInfo','LG_QualityInfo','DataOut','ArousalDat','fitQual','SleepData');
    catch me
        displaytext=me.message;
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
end


displaytext='Complete';
disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

%% Window Selection
function [AnalysisIndex, LGplusinfo,EventsInfo,LG_QualityInfo,DataOut,ArousalDat,fitQual,SleepData] = WindowSelectAndRun(DataEventHypnog_Mat,secslide,VraOn,WindowDuration,ColumnHeads,saveplots)
global settings handletext n winNum
%globals: Fs LGfromFlowVersion rerunspecificwindows plotfigure locttot ignoreCPAPdata;

% *************************************************************************
% ACZ Column Heads
% *************************************************************************
% ColumnHeads=[1    2          4    8      9       10          11       12     13       3        5    6    7];
%   %         [Time RIP_Thorax Flow Hypnog Central Obstructive Hypopnea Desats Arousals RIP_Abdo SpO2 C3A2 Position ]
% *************************************************************************
% *************************************************************************

settings.locttot=0;
settings.ignoreCPAPdata=1;

respwav=DataEventHypnog_Mat(:,ColumnHeads(3));
hypnog=DataEventHypnog_Mat(:,ColumnHeads(4));

if settings.ignoreCPAPdata
    cpap=0*DataEventHypnog_Mat(:,1); %default assume CPAP is absent.
else
    cpap=DataEventHypnog_Mat(:,ColumnHeads(14));
end

position=DataEventHypnog_Mat(:,ColumnHeads(13));
duration=length(respwav);

numwind=floor((duration-WindowDuration*settings.Fs)/(secslide*settings.Fs)); % calculates the number window (of given size, and of given slide) required to encapsulate range
disp(['N possible windows: ' num2str(numwind)])
clear allsleep CPAPoff
if isempty(settings.rerunspecificwindows)
    winnumrange=0:1:numwind-1;
else
    winnumrange=settings.rerunspecificwindows;
end

SleepData=NaN*zeros(numwind,7);
for winNum=winnumrange
    
    % Extract the relevant variables from the particular window.
    hypnogwind=hypnog(winNum*secslide*settings.Fs+1:winNum*settings.Fs*secslide+WindowDuration*settings.Fs);
    CPAPwind=cpap(winNum*secslide*settings.Fs+1:winNum*settings.Fs*secslide+WindowDuration*settings.Fs);
    
    %findest longest durations of wakefulness       
        clear I I2 I3 I4
        I=hypnogwind==4; I2=[I(1)==1;diff(I)]; I3=find(I2==1); I4=find(I2==-1); 
        if length(I4)<length(I3), I4(length(I3),1)=length(hypnogwind)+1; end
        lengthwake=(I4-I3)/settings.Fs;
        if isempty(lengthwake),lengthwake=0; end
        LongestWake=max(lengthwake);
    
    % If it is N1, N2 or N3 sleep the whole window, calculate loop gain:
    allsleep(winNum+1)=(sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0);
    CPAPoff(winNum+1)=max(abs(CPAPwind))<0.5;
    
    FNREM1=sum(hypnogwind==2)/length(hypnogwind);
    FNREM2=sum(hypnogwind==1)/length(hypnogwind);
    FNREM3=sum(hypnogwind==0)/length(hypnogwind);
    FREM=sum(hypnogwind==3)/length(hypnogwind);
    FNREM=FNREM1+FNREM2+FNREM3;
    FWAKE=sum(hypnogwind==4)/length(hypnogwind);
    SleepData(winNum+1,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
    ContainsUnknownSleep=sum(hypnogwind==-1)>0;
    REMdetected=(sum(hypnogwind==3)>0);
    
    if settings.AnalyzeNREMonly
        criteria = (~ContainsUnknownSleep)&(LongestWake<=settings.longestwakeduration)&(REMdetected==0);
    else
        criteria = 1;
    end
    
    %criteria=(sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0)&&(max(abs(CPAPwind))<0.5);
    
    if criteria
        %disp(['Analyzing ' num2str(n) ':' num2str(winNum) '/' num2str(numwind) ', All sleep = ' num2str(allsleep(winNum+1)) ', CPAP off = ' num2str(CPAPoff(winNum+1))]);
        displaytext=['Analyzing ' num2str(n) ':' num2str(winNum) '/' num2str(numwind)];
        disp(displaytext);
        set(handletext,'String',displaytext);
        drawnow;
        AnalysisIndex(winNum+1,:)=[winNum*secslide*settings.Fs+1 winNum*settings.Fs*secslide+WindowDuration*settings.Fs];
        try
            eval(...
                ['[LGplusinfo(winNum+1,:),EventsInfo(winNum+1,:),LG_QualityInfo(winNum+1,:),DataOut{winNum+1},ArousalDat{winNum+1},fitQual{winNum+1}]=' ...
                settings.LGfromFlowVersion '(DataEventHypnog_Mat(winNum*secslide*settings.Fs+1:winNum*settings.Fs*secslide+WindowDuration*settings.Fs,:),ColumnHeads,VraOn,WindowDuration);' ...
                ]);
            %[temp1,temp2,temp3,temp4,temp5,temp6]=LGfromFlowBeta(DataEventHypnog_Mat(winNum*secslide*settings.Fs+1:winNum*settings.Fs*secslide+WindowDuration*settings.Fs,:),ColumnHeads,VraOn,WindowDuration);
        catch me
            disp(['error evaluating LGfromFlow or saving its data: ' me.message])
        end
        
    else
        %disp(['ignore window ' num2str(n) '/' num2str(winNum) ', allsleep=' num2str(allsleep(winNum+1)) ', CPAPoff=' num2str(CPAPoff(winNum+1))]);
        AnalysisIndex(winNum+1,1:2)=NaN;
        LGplusinfo(winNum+1,1:17)=NaN;
        EventsInfo(winNum+1,1:10)=NaN;
        LG_QualityInfo(winNum+1,1:13)=NaN;
        DataOut{winNum+1}=NaN;
        ArousalDat{winNum+1}=NaN;
        fitQual{winNum+1}=NaN;
    end
end
end

