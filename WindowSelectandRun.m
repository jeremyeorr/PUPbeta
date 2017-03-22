function [AnalysisIndex, LGplusinfo,EventsInfo,LG_QualityInfo,DataOut,ArousalDat,fitQual] = WindowTest_LG_Phenotype_publishCode(DataEventHypnog_Mat,windsizeVect,secslide,VraOn,WindowDuration,ColumnHeads,saveplots)
global n winNum LGfromFlowVersion rerunspecificwindows plotfigure;

% *************************************************************************
% ACZ Column Heads
% *************************************************************************
% ColumnHeads=[1    2          4    8      9       10          11       12     13       3        5    6    7];
%   %         [Time RIP_Thorax Flow Hypnog Central Obstructive Hypopnea Desats Arousals RIP_Abdo SpO2 C3A2 Position ]
% *************************************************************************
% *************************************************************************

locttot=0;
ignoreCPAPdata=1;

numWinsize=length(windsizeVect);

respwav=DataEventHypnog_Mat(:,ColumnHeads(3));
hypnog=DataEventHypnog_Mat(:,ColumnHeads(4));

if ignoreCPAPdata
    cpap=0*DataEventHypnog_Mat(:,1); %default assume CPAP is absent.
else
    cpap=DataEventHypnog_Mat(:,ColumnHeads(14));
end

position=DataEventHypnog_Mat(:,ColumnHeads(13));
duration=length(respwav);




for wsv=1:1:numWinsize
    numwind=floor((duration-windsizeVect(wsv)*100)/(secslide*100)); % calculates the number window (of given size, and of given slide) required to encapsulate range
    numwind
    clear allsleep CPAPoff
    if isempty(rerunspecificwindows)
        winnumrange=0:1:numwind-1;
    else
        winnumrange=rerunspecificwindows;
    end
    for winNum=winnumrange
        %for winNum=184:184
        %for winNum=793:1:794
        %for i=0:1:6
        winNum
        % Extract the relevant variables from the particular window.
        hypnogwind=hypnog(winNum*secslide*100+1:winNum*100*secslide+windsizeVect(wsv)*100);
        CPAPwind=cpap(winNum*secslide*100+1:winNum*100*secslide+windsizeVect(wsv)*100);
        position_wind=position(winNum*secslide*100+1:winNum*100*secslide+windsizeVect(wsv)*100);
        % Repair any clipping at the positive or negative rails of the
        % signal
        
        
        
        % If it is N1, N2 or N3 sleep the whole window, calculate Loop
        % Gain:
        allsleep(winNum+1)=(sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0);
        CPAPoff(winNum+1)=max(abs(CPAPwind))<0.5;
        if (sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0)&&(max(abs(CPAPwind))<0.5)% &&(max(abs(CPAPwind-6))<0.5) - (this is accounted for in results)    %&(sum(eventwind)==0) % accounted for in results analysis
            % Convert flow to ventilation signal:
            disp(['analyze window ' num2str(n) '/' num2str(winNum) ', allsleep=' num2str(allsleep(winNum+1)) ', CPAPoff=' num2str(CPAPoff(winNum+1))]);
            AnalysisIndex(winNum+1,:)=[winNum*secslide*100+1 winNum*100*secslide+windsizeVect(wsv)*100];
            try
                eval(...
                    ['[LGplusinfo(winNum+1,:),EventsInfo(winNum+1,:),LG_QualityInfo(winNum+1,:),DataOut{winNum+1},ArousalDat{winNum+1},fitQual{winNum+1}]=' ...
                    LGfromFlowVersion '(DataEventHypnog_Mat(winNum*secslide*100+1:winNum*100*secslide+windsizeVect(wsv)*100,:), ColumnHeads,plotfigure,locttot,VraOn,WindowDuration)' ...
                    ]);
            catch me
                disp(['error evaluating LGfromFlow or saving its data: ' me.message])
            end
            
        else
            disp(['ignore window ' num2str(n) '/' num2str(winNum) ', allsleep=' num2str(allsleep(winNum+1)) ', CPAPoff=' num2str(CPAPoff(winNum+1))]);
            AnalysisIndex(winNum+1,1:2)=NaN;
            LGplusinfo(winNum+1,1:17)=NaN;
            %LGplusinfo_RIP_t(i+1,1:16)=NaN;
            RQAouts(winNum+1,1:11)=NaN;
            RQAouts_Surr(winNum+1,1:11)=NaN;
            Stat(winNum+1,1:5)=NaN;
            Stat_PaO2(winNum+1,1:5)=NaN;
            DQI_0(winNum+1,1:6)=NaN;
            DQI_10(winNum+1,1:6)=NaN;
            EventsInfo(winNum+1,1:10)=NaN;
            EventHistOut{winNum+1}=NaN;
            LG_QualityInfo(winNum+1,1:8)=NaN;
            DataOut{winNum+1}=NaN;
            ArousalDat{winNum+1}=NaN;
            fitQual{winNum+1}=NaN;
        end
    end
end
end




% W_AR1out=[];
% W_AR2out=[];
% W_AR3out=[];
% numWinsize=length(windsizeVect);
%
% respwav=DataEventHypnog_Mat(:,ColumnHeads(3));
% hypnog=DataEventHypnog_Mat(:,ColumnHeads(4));
% duration=length(respwav);
%
% % Fs=100;
% % downSampleFactor=10;
% % Fs=Fs/downSampleFactor;
%
%
% for wsv=1:1:numWinsize
%     numwind=floor((duration-windsizeVect(wsv)*100)/(secslide*100)); % calculates the number window (of given size, and of given slide) required to encapsulate range
%     numwind
%
%     for winNum=0:1:numwind-1
%     %for winNum=793:1:794
%     %for i=0:1:6
%         winNum
%         % Extract the relevant variables from the particular window.
%         hypnogwind=hypnog(winNum*secslide*100+1:winNum*100*secslide+windsizeVect(wsv)*100);
%         % Repair any clipping at the positive or negative rails of the
%         % signal
%
%
%
%         % If it is N1, N2 or N3 sleep the whole window, calculate Loop
%         % Gain:
%         if (sum((hypnogwind~=0)&(hypnogwind~=1)&(hypnogwind~=2))==0)    %&(sum(eventwind)==0) % Edit comment to exclude/include periods with events
%             % Convert flow to ventilation signal:
%             AnalysisIndex(winNum+1,:)=[winNum*secslide*100+1 winNum*100*secslide+windsizeVect(wsv)*100];
%             [LGplusinfo(winNum+1,:),EventsInfo(winNum+1,:),LG_QualityInfo(winNum+1,:),DataOut{winNum+1},ArousalDat{winNum+1}]...
%             = LGfromFlow_plot_V11_0( DataEventHypnog_Mat(winNum*secslide*100+1:winNum*100*secslide+windsizeVect(wsv)*100,:), ColumnHeads, 0 ,0);
%         else
%             AnalysisIndex(winNum+1,1:2)=NaN;
%             LGplusinfo(winNum+1,1:17)=NaN;
%             %LGplusinfo_RIP_t(i+1,1:16)=NaN;
%             RQAouts(winNum+1,1:11)=NaN;
%             RQAouts_Surr(winNum+1,1:11)=NaN;
%             Stat(winNum+1,1:5)=NaN;
%             Stat_PaO2(winNum+1,1:5)=NaN;
%             DQI_0(winNum+1,1:6)=NaN;
%             DQI_10(winNum+1,1:6)=NaN;
%             EventsInfo(winNum+1,1:10)=NaN;
%             EventHistOut{winNum+1}=NaN;
%             LG_QualityInfo(winNum+1,1:8)=NaN;
%             DataOut{winNum+1}=NaN;
%             ArousalDat{winNum+1}=NaN;
%         end
%     end
% end
% end


