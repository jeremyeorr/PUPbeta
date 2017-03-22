function RatingsPUP()
%% Judge quality (0-5 stars)

%% Load data
global AnalyzeDataSpreadsheet handletext settings

if isempty(AnalyzeDataSpreadsheet)
    AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
end

[num,patients] = xlsread(AnalyzeDataSpreadsheet,1,'B3:F5003');
analyzelist = num(:,3);
displaytext='Starting up ratings analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;

%read in savename for ratings
[~,txt,raw] = xlsread(AnalyzeDataSpreadsheet,2,'B3:C20');
settings.savename = char(txt(1,2));
settings.OutputDataDirectory = char(txt(18,2));

if exist([settings.OutputDataDirectory '\' settings.savename '.mat'],'file')==2
    load([settings.OutputDataDirectory '\' settings.savename]);
end

%%
updateratings=0;

for n=1:size(patients,1);
    if analyzelist(n)==0
        displaytext=['Skipping: n=' num2str(n) ', ' char(patients(n,1))];
        disp(displaytext); %set(handletext,'String',displaytext); drawnow;
        continue
    end
    
    try
    LGdata=LGplusinfo{n}(:,7);
    if updateratings&&~(n>size(Ratings,1))&&~isempty(Ratings{n})
        updateratings2=1;
    else
        updateratings2=0;
    end
    
    if updateratings2
        Rating=Ratings{n};
    else
        Rating=NaN;
    end
    
    clear winNum_lastrating
    for winNum=1:length(LGdata)
%         if winNum == length(LGdata)-2
%             kkl=5
%         end
        if ~updateratings2
            Rating(winNum)=NaN;
        end
        figurename=[settings.savename ', n=' num2str(n) ', w=' num2str(winNum-1)]; %added the minus 1 recently
        if ~exist([settings.OutputDataDirectory '\' figurename '.fig'],'file')
            Rating(winNum)=-1;
        end
        if exist([settings.OutputDataDirectory '\' figurename '.fig'],'file')
            try
                if ~isnan(LGdata(winNum,1)) %assumes if col 1 is NaN then all cols are NaN
                    close all force
                    uiopen([settings.OutputDataDirectory '\' figurename '.fig'],1)
                    set(gcf,'units','normalized','outerposition',[0 0 1 1])
                    
                    Rating(winNum)=input('star rating 0-5: ');
%                     Rating(winNum)=7;
                    if exist('winNum_lastrating')&&abs(Rating(winNum)-floor(Rating(winNum)))>0.8&&abs(Rating(winNum)-floor(Rating(winNum)))<1 %the undo button
                        disp('Redo previous, after rating current window')
                        Rating(winNum)=input('star rating 0-5 of current window: ');
                        figurename=[settings.savename ', n=' num2str(n) ', w=' num2str(winNum_lastrating-1)];
                        close all force
                        uiopen([settings.OutputDataDirectory '\' figurename '.fig'],1)
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        Rating(winNum_lastrating)=input('star rating 0-5: ');
                    end
                    winNum_lastrating=winNum;
                end
            catch me
                %Rating(winNum)=NaN;
            end
        end
    end
    Ratings{n}=Rating;
    save([settings.OutputDataDirectory '\' settings.savename],'Ratings','-append');
    cell0='AF12'; %for filename
    xlswrite([settings.OutputDataDirectory '\' char(patients(n,1)) '_results.xlsx'],Ratings{n}',1,cell0);
    catch me
    end
end

%Ratings
%0=no signal 
%1=useless signal
%2=poor signal quality or unusable data for other clear reason e.g. %non-stationarity, mouth breathing, wake/movement, position change
%3=signal is fine but the model fit is not satisfactory / model fit does not capture the essence of the pattern. You would not like to include this data.
%4=satisfactory / model captures the essence of the ventilatory pattern, but not a great fit to the data
%5=excellent / near perfect fit to the data / highly reliable.
%6=example figure quality







