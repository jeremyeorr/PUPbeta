% This version uses LG calculation software neatened for paper publication....
function [LoopGainParameters,EventsInfo,LG_QualityInfo,DataOut,ArousalDat,FitQual] = LGfromFlow_publish(DataEventHypnog_Mat,ColumnHeads,VraSelect,WindowDuration)

global n winNum savename invertflowlist saveplots Pnasaldownisinsp sqrt_scaling usescoredcentralapneas eventsarebreathsfullywithinmargins havescoredcentralhypops
global manualscoringtouchups maxdelaybreaths plotfigure

if ~exist('saveplots','var')
    saveplots=1;
end
if exist('invertflowlist','var')
    if sum(n==invertflowlist)
        invertflow=1;
    else
        invertflow=0;
    end
else
    invertflow=0;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Column heads vector allows backwards compatibility as new channels are
% added to the data hypnog matrix. Channels and labels shown here for
% convenience.
%ColumnHeads=[1      2            4      8        9         10            11         12       13         3           5       6       7];
%            [1=Time 2=RIP_Thorax 3=Flow 4=Hypnog 5=Central 6=Obstructive 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if plotfigure==1
    if 0
        gplot=figure;
    else
        gplot=2;
    end
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Control parameters for operation of algorithm
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if VraSelect==0
    VraOn=0;
else
    VraOn=1;
end
%sqrt_scaling=1;
% mark_onset_of_arousal_only=0;
delete_arousal_in_events=1;
seteventsbasedonvlessthaneupnea2=7;
eventslessthanX=0.7; %default<0.7
findcentralhypopneasandapneas=1; %Turned off on 2013-02-06
removecentralhypops=1;
order=1;
MaxFunEvals=50; %lowered from 500

enforce_eupnea_as_mean=1;
ArousalShift=VraSelect;
polyfitorder=1;
numArousalBreaths=2;

findstartingpoint=0;    %%SS changed to zero on 2013-02-02 since now filtering noise.
detrend_filtfilt_VI=0;  %%SS changed to zero on 2013-02-06 for anatomy measures.

filter_timeconstant=450; %unused.

dt=0.01;
filter_LFcutoff_butter1 = 0.1;
filter_HFcutoff_butter1 = 1;
filter_order1 = 4;
[B_butter1,A_butter1] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Extraction and Preprocessing of flow data for calculation of LG:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Flow=((-1)^Pnasaldownisinsp)*((-1)^invertflow)*DataEventHypnog_Mat(:,ColumnHeads(3)); %%% added invert flow
positionfactor=1; %added 2014-Feb by SS
position_mode=mode(round(positionfactor*DataEventHypnog_Mat(:,ColumnHeads(13)))); %%% added round position
%position_mode=mode(round(DataEventHypnog_Mat(:,ColumnHeads(13)))); %%% added round position
N_position_changes=sum(abs(diff(DataEventHypnog_Mat(:,ColumnHeads(13)))));
Percent_position=mean(position_mode==DataEventHypnog_Mat(:,ColumnHeads(13))); %does this work(?)
Time=DataEventHypnog_Mat(:,ColumnHeads(1));
Arousal=DataEventHypnog_Mat(:,ColumnHeads(9));
hypnog=DataEventHypnog_Mat(:,ColumnHeads(4));
duration=length(Flow);

% Load in the manually scored events.
%obstructive if obstrucive apnea / hypopnea / or mixed apnea

%if using separate channel 14 for mixed events; nb ch14 has been used for CPAP level previously
%Events=DataEventHypnog_Mat(:,ColumnHeads(6))|DataEventHypnog_Mat(:,ColumnHeads(7))|DataEventHypnog_Mat(:,ColumnHeads(14));
%else
Events=DataEventHypnog_Mat(:,ColumnHeads(6))|DataEventHypnog_Mat(:,ColumnHeads(7));

EventsScoredCentralApnea=DataEventHypnog_Mat(:,ColumnHeads(5));
if havescoredcentralhypops
    try
        EventsScoredCentralHyp=DataEventHypnog_Mat(:,ColumnHeads(15));
    catch tempme
        EventsScoredCentralHyp=0*Events;    %row of zeros if no data available
    end
else
    EventsScoredCentralHyp=0*Events;    %row of zeros if no data available
end

[Flow,clipped] = FixClipping_V14_0(Time,Flow,0.05);
if length(Flow)==1
    LoopGainParameters(1:17)=NaN;
    EventsInfo(1:10)=NaN;
    LG_QualityInfo(1:8)=NaN;
    DataOut=NaN;
    ArousalDat=NaN;
    FitQual=NaN;
    return
end

    [time,Vflow,index,time1,VI1,BB_i_start,BB_i_mid,BB_i_end,BB_t,VTi,VTe,VI,VE,Ttot_B,Ti,Te,FLindex1,PeakInspFlow,PeakInspFlow_t,MidInspFlow] = VEfromFlow_sqrt_V14_1(Time,Flow,sqrt_scaling);


if length(VI)==1
    LoopGainParameters(1:17)=NaN;
    %LGplusinfo_RIP_t(1:10)=NaN;
    RQAouts(1:11)=NaN;
    RQAouts_Surr(1:11)=NaN;
    Stat(1:5)=NaN;
    Stat_PaO2(1:5)=NaN;
    DQI_0(1:6)=NaN;
    DQI_10(1:6)=NaN;
    clippedinsp_B=NaN;
    EventsInfo(1:10)=NaN;
    EventHistOut{1}=NaN;
    EventHistOut{2}=NaN;
    LG_QualityInfo(1:8)=NaN;
    DataOut=NaN;
    ArousalDat=NaN;
    FitQual=NaN;
    return
end

%Convert events to breath domain
clippedinsp_B=zeros(length(BB_i_start),1);
clippedexp_B=zeros(length(BB_i_start),1);
AR=zeros(1,length(BB_i_start));
ARonset=zeros(1,length(BB_i_start));
E=zeros(1,length(BB_i_start)); %Events
E_ScoredCentralApnea=zeros(1,length(BB_i_start)); %EventsScoredCentralApnea
E_ScoredCentralHyp=zeros(1,length(BB_i_start)); %EventsScoredCentralHyp
for k=1:length(BB_i_start)
    clippedinsp_B(k)=mean(clipped(BB_i_start(k):BB_i_mid(k)));
    clippedexp_B(k)=mean(clipped(BB_i_mid(k):BB_i_end(k)));
    if 1
        AR(k)=max(Arousal(BB_i_start(k):BB_i_end(k))); %if arousal is within breath
        %E(k)=1-(max(Events(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        if eventsarebreathsfullywithinmargins
            E(k)=1-(min(Events(BB_i_start(k):BB_i_end(k)))); %if event is within breath
            E_ScoredCentralApnea(k)=(min(EventsScoredCentralApnea(BB_i_start(k):BB_i_end(k)))); %if event is within breath
            E_ScoredCentralHyp(k)=(min(EventsScoredCentralHyp(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        else
            E(k)=1-(max(Events(BB_i_start(k):BB_i_end(k)))); %if event is within breath
            E_ScoredCentralApnea(k)=(max(EventsScoredCentralApnea(BB_i_start(k):BB_i_end(k)))); %if event is within breath
            E_ScoredCentralHyp(k)=(max(EventsScoredCentralHyp(BB_i_start(k):BB_i_end(k)))); %if event is within breath
        end
    else
        AR(k)=Arousal(BB_i_start(k)); %if start of breath is in arousal
        E(k)=1-(round(mean(Events(BB_i_start(k):BB_i_end(k))))); %if most of breath is within event
    end
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% Derive Ventilation from thorax and abdomen RIP bands:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
RIP_Thorax=DataEventHypnog_Mat(:,ColumnHeads(2));
RIP_Abdo=DataEventHypnog_Mat(:,ColumnHeads(10));
RIP_Thorax = filtfilt(B_butter1,A_butter1,RIP_Thorax);
RIP_Abdo = filtfilt(B_butter1,A_butter1,RIP_Abdo);


%check for a single shitty RIP trace
RIPinfo='';
if abs(corr(RIP_Thorax,RIP_Abdo))<0.3
    %there is a shitty trace; duplicate the good one.
    if abs(corr(RIP_Abdo,Flow))>abs(corr(RIP_Thorax,Flow));
        RIP_Thorax=RIP_Abdo;
        RIPinfo='Poor Th';
    else
        RIP_Abdo=RIP_Thorax;
        RIPinfo='Poor Ab';
    end
end

VEThorax=zeros(1,length(BB_i_start));
VEAbdo=zeros(1,length(BB_i_start));
for k=1:length(BB_i_start)
    VEThorax(k)=max(RIP_Thorax((BB_i_start(k):BB_i_end(k))))-min(RIP_Thorax((BB_i_start(k):BB_i_end(k))));
    VEAbdo(k)=max(RIP_Abdo((BB_i_start(k):BB_i_end(k))))-min(RIP_Abdo((BB_i_start(k):BB_i_end(k))));
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Take ventilation as the average of inspired/expired levels, and find
% RIP thorax to abdo ratios
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%Take ventilation as the average of inspired/expired levels.
VI=(VI+VE)/2; % Mean of inspiration and expiration volumes
VI=VI/mean(VI); %normalize
VEThorax=VEThorax/mean(VEThorax); %normalize
VEAbdo=VEAbdo/mean(VEAbdo); %normalize
VEThorax(VEThorax>3)=3;
VEAbdo(VEAbdo>3)=3;

%THABratio
alphas=[0:0.01:1];
SSE_alphas=zeros(1,length(alphas));
corr_alphas=zeros(1,length(alphas));
for k=1:length(alphas)
    RIP=alphas(k)*VEThorax+(1-alphas(k))*VEAbdo;
    m=VI\RIP;
    Yexpected=RIP*m;
    SSE_alphas(k)=sum((VI-Yexpected).^2);
    corr_alphas(k)=corr(VI',RIP');
end
[~,tempi]=min(SSE_alphas);
RIPalphaSSE=alphas(tempi);
[~,tempi]=max(corr_alphas);
RIPalphaCorr=alphas(tempi);

VI_RIP1=RIPalphaCorr*VEThorax+(1-RIPalphaCorr)*VEAbdo;
VI_RIP2=RIPalphaSSE*VEThorax+(1-RIPalphaSSE)*VEAbdo;
m1=VI'\VI_RIP1'; m2=VI'\VI_RIP2';

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calculate LG from ventilation data:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
TimeB = BB_t-BB_t(1); % The time of each breath, from the first breath....[breath time series]

% *********************************************************************
% if there is a missing section of samples either at the end or the
% beginning of the period, break.
% *********************************************************************
if TimeB(end)<0.8*WindowDuration
    LoopGainParameters(1:17)=NaN;
    %LGplusinfo_RIP_t(1:16)=NaN;
    RQAouts(1:11)=NaN;
    RQAouts_Surr(1:11)=NaN;
    Stat(1:5)=NaN;
    Stat_PaO2(1:5)=NaN;
    DQI_0(1:6)=NaN;
    DQI_10(1:6)=NaN;
    clippedinsp_B=NaN;
    EventsInfo(1:10)=NaN;
    EventHistOut{1}=NaN;
    EventHistOut{2}=NaN;
    LG_QualityInfo(1:8)=NaN;
    DataOut=NaN;
    ArousalDat=NaN;
    FitQual=NaN;
    return
end
% *********************************************************************

% *********************************************************************
% Detrend the ventilation signal using a FIRST order polynomial fit
% *********************************************************************
if detrend_filtfilt_VI
    polyfitresults=polyfit(DataEventHypnog_Mat(BB_i_start,1)',VI,polyfitorder);
    VImovingeupnea=0*VI;
    for k=1:(polyfitorder+1)
        VImovingeupnea=VImovingeupnea+polyfitresults(k)*DataEventHypnog_Mat(BB_i_start,1)'.^(polyfitorder+1-k);
    end
    if 0
        polyfitresultsVth=polyfit(DataEventHypnog_Mat(BB_i_start,1)',VEThorax,polyfitorder);
        VImovingeupneaTH=0*VEThorax;
        for k=1:(polyfitorder+1)
            VImovingeupneaTH=VImovingeupneaTH+polyfitresultsVth(k)*DataEventHypnog_Mat(BB_i_start,1)'.^(polyfitorder+1-k);
        end
        polyfitresultsVab=polyfit(DataEventHypnog_Mat(BB_i_start,1)',VEAbdo,polyfitorder);
        VImovingeupneaAB=0*VEAbdo;
        for k=1:(polyfitorder+1)
            VImovingeupneaAB=VImovingeupneaAB+polyfitresultsVab(k)*DataEventHypnog_Mat(BB_i_start,1)'.^(polyfitorder+1-k);
        end
        VEThorax=VEThorax+mean(VEThorax)-VImovingeupneaTH;
        VEAbdo=VEAbdo+mean(VEAbdo)-VImovingeupneaAB;
    end
    
    VI=VI+mean(VI)-VImovingeupnea;
end
%re-normalize:
VI=VI/mean(VI); Veupnea=1;
VEThorax=VEThorax/mean(VEThorax); %normalize thorax RIP
VEAbdo=VEAbdo/mean(VEAbdo); %normalize abdo RIP
% *********************************************************************

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Update events based on eupnea, and detection of "detected" central
% apneas... i.e. any breath greater than eupnoea has events cleared,
% any breath at less than 70% of eupnoea is set as an event. Also, any
% breath detected as a central event (depending on criteria) has event
% cleared for purpose of model fitting.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    plot_on2=plotfigure;
else
    plot_on2=0;%Override - don't plot plots from ObstructBreathScore_V14_0
end

[E1,Ecentralapnea,Ecentralhypop] = ObstructBreathScore_V14_1(Veupnea, E, VI, VEThorax, VEAbdo, BB_i_start, removecentralhypops, findcentralhypopneasandapneas,seteventsbasedonvlessthaneupnea2,eventslessthanX, plot_on2);
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if usescoredcentralapneas %Use scored central events (if they exist)
    E1(E_ScoredCentralHyp|E_ScoredCentralApnea)=1; %1=unobstructed
    Ecentralapnea(E_ScoredCentralApnea==1)=1; %1=central apnea
    Ecentralhypop(E_ScoredCentralHyp==1)=1; %1=central hypopnea
end

if 0
    %try additional events manipulation based on flow trace:
    %reset:
    E1(VI<0.7)=0;
    E1(VI>0.7)=1;
    E1(Ecentralapnea==1)=1;
    E1(Ecentralhypop==1)=1;
    %override based on FLindex and Ti/Ttot
    c_=0.3; b_=0.5;
    E1(VI>1&VI<((1-c_)/b_*FLindex1+c_)&Ti/Ttot_B>0.5)=0; %event %check greater/less thans...
    E1(VI<1&VI<((1-c_)/b_*FLindex1+c_))=0; %event
    E1(VI>1.3)=1; %not event
    E1(VI<1.1&FLindex1>0.05&Ti/Ttot_B>0.6)=0;
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Alternative events, purely for the purpose of arousal marking
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
E_recover=zeros(1,length(VI));
E_Terminate=zeros(1,length(VI));
E_recover(VI>1)=1;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if 0
    if plotfigure==1
        figure(100),
        ax100(1)=subplot(311);plot(DataEventHypnog_Mat(BB_i_start,1),VI,DataEventHypnog_Mat(BB_i_start,1),[zeros(1,length(BB_i_start));ones(1,length(BB_i_start))],':k');
        ax100(2)=subplot(312);plot(DataEventHypnog_Mat(BB_i_start,1),VEThorax,DataEventHypnog_Mat(BB_i_start,1),VEAbdo,DataEventHypnog_Mat(BB_i_start,1),[zeros(1,length(BB_i_start));ones(1,length(BB_i_start))],':k');
        ax100(3)=subplot(313);plot(DataEventHypnog_Mat(BB_i_start,1),Ecentralapnea,'.',DataEventHypnog_Mat(BB_i_start,1),Ecentralhypop,'.');
        linkaxes(ax100,'x');
    end
end


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This section identifies the duration of each arousal and links it to
% the breath domain arousal markings
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR_Rise=AR(2:end)-AR(1:end-1);
AR_Rise(AR_Rise<1)=0;
ArousalDuration=zeros(1,length(AR));
for jk=1:1:length(AR_Rise)
    if AR_Rise(jk)==1
        time_start=BB_i_start(jk+1);
        time_count=time_start;
        while Arousal(time_count)==0&&time_count<length(Time)
            time_count=time_count+1;
        end
        while Arousal(time_count)==1&&time_count<length(Time)
            time_count=time_count+1;
        end
        ArousalDuration(jk+1)=Time(time_count)-Time(time_start);
        
        %AR(k)=max(Arousal(BB_i_start(k):BB_i_mid(k)));
    end
end

ArousalCount=1;
for jk=2:1:length(Arousal)
    if Arousal(jk)==1&&Arousal(jk-1)==0 % Arousal Start
        ArousalDat(ArousalCount,1)=Time(jk);
    elseif Arousal(jk)==0&&Arousal(jk-1)==1 %Arousal End
        ArousalDat(ArousalCount,2)=Time(jk);
        ArousalCount=ArousalCount+1;
    end
end
if ArousalCount==1
    ArousalDat(ArousalCount,2)=NaN;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%This section places notations on event markings to allow
%implementation of the desired VRA modelling and data fitting.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR1=AR;
if delete_arousal_in_events
    [ AR1, numArousalBreaths, E_Terminate] = ArousalShift_V14_0b(ArousalShift,AR1,E1,E_Terminate,BB_i_start,numArousalBreaths);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Run the autoregressive loop gain calculation for a range of delays,
% and pick the one that results in the best (lowest) cost function.
%**********************************************************************
%**********************************************************************
forcedelays=[1:maxdelaybreaths NaN]; % delay from 1 to 5 breaths.
Nforcedelays=length(forcedelays)-1; % number of iterations of delay before searching for the lowest resultant cost function.


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Identify the best starting breath from the window, based on: starting
% after the maximum possible delay, and being the start of the longest
% section of data in the first 2 minutes that is costed in the AR
% function (i.e: no event)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Ttot_ser=(diff(TimeB)); % T-tot series
Ttot_med=median(Ttot_ser);
Ttot_ser_prev=[Ttot_med Ttot_ser(1:end)];

while 1
    Ventilation_EventScoring_data=[VI'-Veupnea E1' AR1(1:length(E1))' TimeB' Ecentralapnea' [Ttot_ser.'; 0] Ttot_ser_prev.'];
    
    % Ttot_ser_prev=[Ttot_med Ttot_ser(1:end-1)];
    % Ventilation_EventScoring_data=[VI'-Veupnea E1' AR1(1:length(E1))' TimeB' Ecentralapnea' [Ttot_ser.'; 0] [Ttot_ser_prev.';0]];
    
    %find arousal onset
    ARonset=AR1;
    ARonset(1:end)=0;
    for i=2:1:length(AR1)
        if AR1(i)==1&&AR1(i-1)==0
            ARonset(i)=1;
        end
    end
    
    %counting events, note this only counts obstructive events.
    N_events=0;
    % duration=0;
    for k=2:length(E)
        if E(k)==0&&E(k-1)==1
            N_events=N_events+1;
            %         duration(N_events)=duration(N_events-1)
        end
    end
    
    if  N_events<1
        LoopGainParameters(1:17)=NaN;
        %LGplusinfo_RIP_t(1:16)=NaN;
        RQAouts(1:11)=NaN;
        RQAouts_Surr(1:11)=NaN;
        Stat(1:5)=NaN;
        Stat_PaO2(1:5)=NaN;
        DQI_0(1:6)=NaN;
        DQI_10(1:6)=NaN;
        clippedinsp_B=NaN;
        EventsInfo(1:10)=NaN;
        EventHistOut{1}=NaN;
        EventHistOut{2}=NaN;
        LG_QualityInfo(1:8)=NaN;
        DataOut=NaN;
        ArousalDat=NaN;
        FitQual=NaN;
        disp('No scored events, E=0');
        return
    end
    
    % Ventilation_EventScoring_data(end,:)=[];
    %ARonset(end,:)=[];
    modelpolyfitorder=3;
    [Error,Vdr_est,VAr_est,LoopGainParameters,Optimal_vent_control_parameters,Optimal_MSE,FitQual,i_1,i_end] = FindBestModelParameters(Ventilation_EventScoring_data,VraOn,Veupnea,modelpolyfitorder);
    
    
    
    
    
    if isnan(LoopGainParameters(1))==1
        LoopGainParameters(1:17)=NaN;
        %LGplusinfo_RIP_t(1:16)=NaN;
        RQAouts(1:11)=NaN;
        RQAouts_Surr(1:11)=NaN;
        Stat(1:5)=NaN;
        Stat_PaO2(1:5)=NaN;
        DQI_0(1:6)=NaN;
        DQI_10(1:6)=NaN;
        clippedinsp_B=NaN;
        EventsInfo(1:10)=NaN;
        EventHistOut{1}=NaN;
        EventHistOut{2}=NaN;
        LG_QualityInfo(1:8)=NaN;
        DataOut=NaN;
        ArousalDat=NaN;
        FitQual=NaN;
        return
    end
    
    Arthressample=Vdr_est+Veupnea;
    Arthressample([1:i_1 i_end:end])=[];
    ARonset([1:i_1 i_end:end])=[];
    Arthressample(ARonset==0)=[];
    Arthres=mean(Arthressample);
    N_arousals_used=sum(ARonset);
    
    
    %% Plot
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %If plotfigure=1, create plot
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (plotfigure==1)
        
        if 0
            plotflow=Flow; plotflow_t=Time;
        else
            plotflow=Vflow; plotflow_t=time;
        end
        if 1
            RIP_Thorax_plot=DataEventHypnog_Mat(:,ColumnHeads(2));
            RIP_Abdo_plot=DataEventHypnog_Mat(:,ColumnHeads(10));
            RIP_Thorax_plot = filtfilt(B_butter1,A_butter1,RIP_Thorax_plot);
            RIP_Abdo_plot = filtfilt(B_butter1,A_butter1,RIP_Abdo_plot);
        end
        set(gcf,'Color',[1 1 1]);
        ax1(1)=subplot(311);
        plot(BB_t(AR1>0),AR1(AR1>0),'r.',BB_t(E1<1),E1(E1<1),'k.',BB_t(E<1),E(E<1)+1,'c.',BB_t(Ecentralhypop>0),Ecentralhypop(Ecentralhypop>0)-2,'b.',...
            BB_t(Ecentralapnea>0),Ecentralapnea(Ecentralapnea>0)-3,'g.');
        set(gca,'YColor',[1 1 1],'XColor',[1 1 1],...
            'FontSize',7,...
            'Position',[0.131 0.802 0.775 0.031],...
            'FontName','Arial');
        xlim([min(Time) max(Time)]);
        title([num2str(n) ' ' num2str(winNum) ' ' RIPinfo],'fontname','Arial Narrow');
        
        ax1(2)=subplot(312);
        alpha_=99.5;
        nflow=1.5*(prctile(plotflow,alpha_)-prctile(plotflow,100-alpha_));    mflow=mean(plotflow);
        nTH=2*(prctile(RIP_Thorax_plot,alpha_)-prctile(RIP_Thorax_plot,100-alpha_));    mTH=mean(RIP_Thorax_plot);
        nAB=2*(prctile(RIP_Abdo_plot,alpha_)-prctile(RIP_Abdo_plot,100-alpha_));     mAB=mean(RIP_Abdo_plot);
        plot(plotflow_t,0.5+(plotflow-mflow)/nflow,'k',...
            Time,-0.25+(RIP_Thorax_plot-mTH)/nTH,'k',...)
            Time,-1+(RIP_Abdo_plot-mAB)/nAB,'k');
        set(gca,'YTick',[-1 -0.25 0.5],...
            'XColor',[1 1 1],...
            'Position',[0.131 0.455 0.775 0.321],'FontName','Arial','FontSize',7);
        box('off');
        ylim([-1.6 1.6]);
        set(gca,'YDir','normal','FontSize',7);
        set(gca,'YTickLabel',{'Abdo RIP', 'Thorax RIP','Flow'});
        box('off');
        xlim([min(Time) max(Time)]);
        hold('on')
        if 1
            %plot EEG and arousals scored over the top in green
            EEGC3A2=DataEventHypnog_Mat(:,ColumnHeads(12))-mean(DataEventHypnog_Mat(:,ColumnHeads(12)));
            plot(Time,1.2+0.2*EEGC3A2/(prctile(EEGC3A2,99.5)-prctile(EEGC3A2,0.05)),'k');
            tempds=4; tempsize=0.2; tempoffset=1.2-0.1;
            tempx=downsample(Time,tempds);
            tempy=tempsize*downsample(Arousal,tempds)+tempoffset;
            tempxx=[tempx;flipud(tempx)];
            tempyy=[tempy;tempoffset+zeros(length(tempy),1)];
            fill(tempxx,tempyy,[0 1 0],'EdgeColor','none','FaceAlpha',0.75);
            %plot final event breaths over the top of flow in blue
            dt_new=0.25;
            startT=ceil(BB_t(1));
            endT=ceil(BB_t(length(BB_t)));
            time_dt=startT:dt_new:endT;
            E1_rs=0*time_dt;
            for i=1:length(time_dt);
                E1_rs(i) = 1-E1(find(Time(BB_i_start)<=time_dt(i),1,'last'));
            end
            C1_rs=0*time_dt;
            for i=1:length(time_dt);
                C1_rs(i) = Ecentralapnea(find(Time(BB_i_start)<=time_dt(i),1,'last'));
            end
            E1_rs=E1_rs&(~C1_rs);
            tempsize=0.25; tempoffset=0.375;
            tempy=tempsize*(E1_rs)+tempoffset;
            tempxx=[time_dt fliplr(time_dt)];
            tempyy=[tempy tempoffset+zeros(1,length(tempy))];
            fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
            hold('on')
            tempsize=0.25; tempoffset=0.375;
            tempy=tempsize*(C1_rs)+tempoffset;
            tempxx=[time_dt fliplr(time_dt)];
            tempyy=[tempy tempoffset+zeros(1,length(tempy))];
            fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
            %figure(980); plot(time_dt,E1_rs,'r',BB_t,E1,'r.')
        end
        
        
        ax1(3)=subplot(313);
        hold('on');
        %plot event breaths over the top of flow in blue
        tempsize=0.5; tempoffset=0.75;
        tempy=tempsize*(E1_rs)+tempoffset;
        tempxx=[time_dt fliplr(time_dt)];
        tempyy=[tempy tempoffset+zeros(1,length(tempy))];
        fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
        tempsize=0.5; tempoffset=0.75; %updated 6/2014
        tempy=tempsize*(C1_rs)+tempoffset;
        tempxx=[time_dt fliplr(time_dt)];
        tempyy=[tempy tempoffset+zeros(1,length(tempy))];
        fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
        
        stairs(BB_t(1:i_end),Ventilation_EventScoring_data(1:i_end,1)+Veupnea,'color',[0.4 0.4 0.7]);
        plot(BB_t(1:i_end),[zeros(1,length(BB_t(1:i_end))); ones(1,length(BB_t(1:i_end)))],'k:');
        stairs(BB_t(i_1:i_end)+0.3,Vdr_est(i_1:end)'+Veupnea+VAr_est(i_1:i_end)','color',[0.1 0.8 0.1]);
        plot(BB_t(i_1:i_end),Vdr_est(i_1:end)'+Veupnea,'color',[0 0 0]);
        %plot(BB_t(i_1:i_end),W(i_1:end)','r');
        %plot(BB_t(i_1:i_end),VAr_est(i_1:i_end),'b');
        %xlabel(num2str(LGplusinfo(4),2))
        %LGplusinfo=[Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA_est];
        xlabel(['RsqTotal=' num2str(FitQual(2)) ' F=' num2str(Optimal_MSE) ' LG1=' num2str(LoopGainParameters(6),2) ' Tn=' num2str(LoopGainParameters(5),2) ' LGn=' num2str(LoopGainParameters(4),2) ' delay=' num2str(LoopGainParameters(3),2) ' ArThr' num2str(Arthres,2) ])
        hold('off');
        set(gca,'FontSize',7,...
            'Position',[0.131 0.221 0.775 0.216],...
            'FontName','Arial');
        box('off');
        xlim([min(Time) max(Time)]);
        linkaxes(ax1,'x');
        if manualscoringtouchups
            disp('click left of data if no touch ups needed, or left then right to add events in range, or right then left to remove events in range');
            [t(1),~,~]=ginput(1);
            if t(1)<min(Time)
                manualscoringtouchupsincomplete=0;
            else
                [t(2),~,button]=ginput(1);
                manualscoringtouchupsincomplete=1;
                if button==1
                if t(1)<t(2)
                    E1(BB_t>t(1)&BB_t<t(2))=0; %add events...
                elseif t(1)>t(2)
                    E1(BB_t>t(2)&BB_t<t(1))=1; %remove events...
                end
                elseif button>1 %right click
                if t(1)<t(2)
                    Ecentralapnea(BB_t>t(1)&BB_t<t(2))=1; %add events...
                elseif t(1)>t(2)
                    Ecentralapnea(BB_t>t(2)&BB_t<t(1))=0; %remove events...
                end  
                end
            end
        else
            manualscoringtouchupsincomplete=0;
        end
    else
        manualscoringtouchupsincomplete=0;
    end
    
    if manualscoringtouchupsincomplete==0
        break
    end
end %end while manualscoringtouchups loop

if saveplots
    saveas(gcf,[savename ', n=' num2str(n) ', w=' num2str(winNum)], 'fig')
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%LGplusinfo=[Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA_est];
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% LG output vector which includes some other key results parameters.
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
LoopGainParameters=[BB_t(i_1) LoopGainParameters(1:2) 0 LoopGainParameters(4:7) LoopGainParameters(3) LoopGainParameters(8) 0 Arthres Optimal_MSE mean(Ttot_B) std(Ttot_B) median(Ttot_B) iqr(Ttot_B)];
LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position 0];
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%'LGplusinfo=[breath Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA1_est VRA2 Arthres1 Optimal_MSE mean(Ttot_B) std(Ttot_B) median(Ttot_B) iqr(Ttot_B)];'
% old
%'LGplusinfo=[breath Gtot_est tau1_est delay_est LG_n_est T_n_est LG_1_est LG_2_est VRA_est Arthres Optimal_MSE mean(Ttot_B) std(Ttot_B) median(Ttot_B) iqr(Ttot_B)]'
% Current.

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calculate the Events Information Parameters:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
E=Ventilation_EventScoring_data(:,2);
AR=Ventilation_EventScoring_data(:,3);
eventHist=zeros(length(E),1);
eventCurrent=0;
eventCount=0;
for b=1:1:length(E)
    if eventCurrent==0
        if E(b)==0
            eventCurrent=1;
            eventCount=eventCount+1;
        elseif E(b)==1
            
        end
    elseif eventCurrent==1
        if E(b)==0
            eventCount=eventCount+1;
        elseif E(b)==1
            eventHist(eventCount)=eventHist(eventCount)+1;
            eventCurrent=0;
            eventCount=0;
        end
    end
    if (b==length(E))&&(eventCurrent==1)
        eventHist(eventCount)=eventHist(eventCount)+1;
    end
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Breath by Breath Events:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR_Rise=AR(i_1+1:end)-AR(i_1:end-1);
AR_Rise(AR_Rise<1)=0;
Eup=1-E;
EventRise=Eup(i_1+1:end)-Eup(i_1:end-1);
EventRise(EventRise<1)=0;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Time Based Events:
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
AR_Time=Arousal(ceil(TimeB(i_1)*dt)+1:ceil(TimeB(i_1)*dt+600*dt))-Arousal(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt-1));
AR_Time(AR_Time<1)=0;
EventTime=Events(ceil(TimeB(i_1)*dt)+1:ceil(TimeB(i_1)*dt+600*dt))-Events(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt-1));
EventTime(EventTime<1)=0;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

EventsInfo=[length(Events)/100 mean(Arousal(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt))) sum(AR_Time)  mean(Events(ceil(TimeB(i_1)*dt):ceil(TimeB(i_1)*dt+600*dt)))  sum(EventTime) length(VI) mean(AR(i_1+1:end)) sum(AR_Rise)  mean(Eup)  sum(EventRise) ];
EventHistOut{1}=eventHist;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
%     size(BB_i_start(i_1:i_end))
%     size(BB_i_mid(i_1:i_end))
%     size(BB_i_end(i_1:i_end))
%     size(E1(i_1:i_end))
%     size(E_recover(i_1:i_end))
%     size(W(i_1:i_end))
%     size(Vdr_est(i_1:i_end))
%
%DataOut=[Time(BB_i_start(i_1:i_end)) Time(BB_i_mid(i_1:i_end)) Time(BB_i_end(i_1:i_end)) VI(i_1:i_end).' Vdr_est(i_1:i_end).' E1(i_1:i_end).' E_recover(i_1:i_end).' E_Terminate(i_1:i_end).' W(i_1:i_end).' ];
DataOut=[Time(BB_i_start(i_1:i_end)) Time(BB_i_mid(i_1:i_end)) Time(BB_i_end(i_1:i_end)) VI(i_1:i_end).' Vdr_est(i_1:i_end) E1(i_1:i_end).' E_recover(i_1:i_end).' E_Terminate(i_1:i_end).' Error(i_1:i_end) ];
%ArousalDat
% bigsize=size(LoopGainParameters)
% bigsize=size(FitQual)
% bigsize=size(EventsInfo)
% bigsize=size(LG_QualityInfo)

end
