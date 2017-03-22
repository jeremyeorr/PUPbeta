function [time,VFlow_filtered1,index,time1,VI1,BB_i_start,BB_i_mid,BB_i_end,BB_t,VTi,VTe,VI,VE,Ttot,Ti,Te,FLindex1,PeakInspFlow,PeakInspFlow_t,MidInspFlow] = VEfromFlow_sqrt_V14_1(timeScale,respwav,sqrt_scaling) 

minimum_figs=1;
Vflow=respwav;
time=timeScale;

N = length(Vflow);
dt=(time(end)-time(1))/(length(time)-1);

% Offset such that the mean nasal pressure signal=0
detrend_flow=1;
if detrend_flow
    Vflow=Vflow-mean(Vflow);
end

% Square root scaling of the presure signal to make it better approximate a
% true flow signal:
if sqrt_scaling
    Vflow(Vflow>0)=Vflow(Vflow>0).^0.5;
    Vflow(Vflow<0)=-((-Vflow(Vflow<0)).^0.5);
end

T0 = time(end)-time(1); %may be redundant


indexstart=0;
indexend=0;
if 1 %clip trace to first and last zero crossings
    for i=2:length(Vflow),
        if (Vflow(i)>0)&&(Vflow(i-1)<=0)
            indexstart=i;
            break
        end
    end
    for i=length(Vflow):-1:2
        if (Vflow(i)>0)&&(Vflow(i-1)<=0)
            indexend=i;
            length(Vflow)-indexend;
            break
        end
    end
    
    if indexend<=indexstart
        % Only one zero crossing - set all outputs to NaN, and return
        % values. 
        time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN;time=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; VFlow_filtered1=NaN; index=NaN;
        return;
    else
        Vflow(indexend:end)=[];
        Vflow(1:indexstart)=[];
        time(indexend:end)=[];
        time(1:indexstart)=[];
        %time=time(indexstart:indexend);
        T0 = time(end)-time(1);
        % if detrend_flow
        %     Vflow=Vflow-mean(Vflow);
        % end
    end
    

end

if length(Vflow)<1000 % At least 10 secs of data
    time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN;time=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; VFlow_filtered1=NaN; index=NaN;
    return;
end

%Filter gently for use in analysis
filter_LFcutoff_butter1 = 0.1;
filter_HFcutoff_butter1 = 3;
filter_order1 = 2;
[B_butter1,A_butter1] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));
VFlow_filtered1 = filtfilt(B_butter1,A_butter1,Vflow);

%% Filter strongly to aid peak identification
[Presp,F1] = pwelch(VFlow_filtered1-mean(VFlow_filtered1),length(Vflow),0,length(Vflow),1/dt);
df=1/T0;
guess_resp_freq=1/4;
Fi=floor((guess_resp_freq/2)/df+1); %25 sec breaths
Fii=ceil((guess_resp_freq*2)/df+1); %0.5 sec breaths
[temp,index_f]=max(Presp(Fi:Fii));
index_f2=index_f+Fi-1;

%peak fit quadratic:
range3 = (index_f2-1):(index_f2+1);
[fresp_peak,Presp_peak] = PeakFitQuadratic_V14_0(F1(range3),Presp(range3));

%Filter strongly to aid peak identification, using LFcutoff as half the expected,
%and HF cutoff as twice the expected
filter_LFcutoff_butter2 = fresp_peak*.2;        
filter_HFcutoff_butter2 = fresp_peak*1.4;
FiltVect=[filter_LFcutoff_butter2 filter_HFcutoff_butter2]/(1/dt/2);
FiltVect(FiltVect<=0)=0.000001;
FiltVect(FiltVect>=1)=0.999999;

filter_order2 = 2;
[B_butter2,A_butter2] = butter(filter_order2,FiltVect);
VFlow_filtered2 = filtfilt(B_butter2,A_butter2,Vflow);

VFlow_filtered_threshold = 0; %-Inf

% if ~minimum_figs   
% figure(200), plot(time,VFlow_filtered1,time,VFlow_filtered2)
% end


if length(Vflow)<10;
    % Only one zero crossing - set all outputs to NaN, and return
    % values.
    time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN;time=NaN; VFlow_filtered1=NaN; index=NaN;
    return;
end
%%detecting max flow values (i.e. expiration..)
clear VFlow_crossing VFlow_crossing_i VFlow_crossing_t;
count=1;
for i=1:length(Vflow),
    if (i>1)&&(i<length(Vflow))
        if (VFlow_filtered2(i)>VFlow_filtered2(i-1))&&(VFlow_filtered2(i)>VFlow_filtered2(i+1))&&(VFlow_filtered2(i)>VFlow_filtered_threshold)
            VFlow_crossing_i(count)=i;
            VFlow_crossing(count)=VFlow_filtered1(VFlow_crossing_i(count));
            VFlow_crossing_t(count)=time(VFlow_crossing_i(count));
            count=count+1;
            
        end
    end
end
if count==1
    time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN;time=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; VFlow_filtered1=NaN; index=NaN;
    return;
end
if 0
%remove last two crossings... always get me in trouble!
VFlow_crossing_i(end)=[];
VFlow_crossing(end)=[];
VFlow_crossing_t(end)=[];
% VFlow_crossing_i(count-2)=[];
% VFlow_crossing(count-2)=[];
% VFlow_crossing_t(count-2)=[];
end

if ~minimum_figs   
figure(3), subplot(3,1,1),plot(time,VFlow_filtered1,time,VFlow_filtered2,VFlow_crossing_t,VFlow_crossing,'.')
end
%Store data
VFlow_crossing_step1_i=VFlow_crossing_i;

% start from peak data and move leftwards to find zero crossing:
for count=1:length(VFlow_crossing_i)
    for k=1:999999
        if (VFlow_filtered1(VFlow_crossing_i(count))>=0)&&(VFlow_filtered1(VFlow_crossing_i(count)-1)<0)
            break;
        else
            if VFlow_crossing_i(count)==2 %hitting left hand limit
                break
            else
                VFlow_crossing_i(count)=VFlow_crossing_i(count)-1;%moving left...
            end
        end
    end
    VFlow_crossing(count)=VFlow_filtered1(VFlow_crossing_i(count));
    VFlow_crossing_t(count)=time(VFlow_crossing_i(count));
end

VFlow_crossing(1)=[];
VFlow_crossing_t(1)=[];
VFlow_crossing_i(1)=[];

i_=0;
for i=1:length(VFlow_crossing)
    if i>1
        if VFlow_crossing_i(i-i_)==VFlow_crossing_i(i-i_-1)
            VFlow_crossing(i-i_)=[];
            VFlow_crossing_t(i-i_)=[];
            VFlow_crossing_i(i-i_)=[];
            i_=i_+1;
        end
    end
end

if ~minimum_figs   
figure(3), subplot(3,1,2),plot(time,Vflow,VFlow_crossing_t,VFlow_crossing,'.r')
end
%% Shift from start inspiration, head rightwards to find end inspiration

VFlow_crossingE=[]; VFlow_crossingE_t=[]; VFlow_crossingE_i=[];
VFlow_crossingE_i = VFlow_crossing_i+2;

for count=1:length(VFlow_crossing_i)
    for k=1:999999
        if (VFlow_filtered1(VFlow_crossingE_i(count))<0)&&(VFlow_filtered1(VFlow_crossingE_i(count)-1)>=0)||(VFlow_crossingE_i(count)>=length(VFlow_filtered1)) %last OR ensures you don't exceed the length of the trace 
            break;
        else
            if VFlow_crossingE_i(count)==2
                break
            else
                VFlow_crossingE_i(count)=VFlow_crossingE_i(count)+1; %increase index
            end
        end
    end
    VFlow_crossingE(count)=VFlow_filtered1(VFlow_crossingE_i(count));
    VFlow_crossingE_t(count)=time(VFlow_crossingE_i(count));
end
%         length(VFlow_crossing_i)
%         length(VFlow_crossingE_i)
%         length(VFlow_crossingE)

if length(VFlow_crossingE)==0
    time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; time=NaN; VFlow_filtered1=NaN; index=NaN;
    return;
end

i_=0;
for i=1:length(VFlow_crossingE)
    if i>1
        % If the crossing index is the same as the previous one, then need
        % to delet the current crossing, and therefore, also the associated
        % inspiratory crossing:
        if VFlow_crossingE_i(i-i_)==VFlow_crossingE_i(i-i_-1)
            %i
            VFlow_crossingE(i-i_)=[];
            VFlow_crossingE_t(i-i_)=[];
            VFlow_crossingE_i(i-i_)=[];
            
            % Double check.......
            % ************************************************************
            VFlow_crossing(i-i_)=[];
            VFlow_crossing_t(i-i_)=[];
            VFlow_crossing_i(i-i_)=[];
            % ************************************************************
            
            i_=i_+1;
        end
    end
end
%         length(VFlow_crossing_i)
%         length(VFlow_crossingE_i)
%         length(VFlow_crossingE)

if ~minimum_figs   
figure(3), subplot(3,1,3),plot(time,VFlow_filtered1,VFlow_crossing_t,VFlow_crossing,'.r',VFlow_crossingE_t,VFlow_crossingE,'.k')
end

%% Separate flow into positive and negative for integration later.

VFlow_filtered1_neg = ones(1,length(VFlow_filtered1)); %initialise for speed
VFlow_filtered1_pos = VFlow_filtered1_neg; %initialise for speed
for i=1:length(VFlow_filtered1)
    if VFlow_filtered1(i)>=0
        VFlow_filtered1_pos(i) = VFlow_filtered1(i);
        VFlow_filtered1_neg(i) = 0;
    else
        VFlow_filtered1_neg(i) = VFlow_filtered1(i);
        VFlow_filtered1_pos(i) = 0;
    end
end

%% First run Ttot data
if 1 %insert breaths when not detected...
    clear Ttot BB_t BB_i Ti Te BB_i_mid BB_i_start BB_i_end VTi VTe;

    if length(VFlow_crossing)<2
        % Only one zero crossing - set all outputs to NaN, and return
        % values.
        time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; time=NaN; VFlow_filtered1=NaN; index=NaN;
        return;
    end
    for i=1:length(VFlow_crossing)
        if i<length(VFlow_crossing)
            %         length(time)
            %         length(VFlow_crossing_i)
            %         length(VFlow_crossingE_i)
            %         length(VFlow_crossing)
            %         i
            Ttot_temp(i)=time(VFlow_crossing_i(i+1))-time(VFlow_crossing_i(i));
            Ti_temp(i)=time(VFlow_crossingE_i(i+1))-time(VFlow_crossing_i(i));
        end
    end
    
    medianTtot1 = median(Ttot_temp);
    medianTi1 = median(Ti_temp);
    medianTe1 = medianTtot1-medianTi1;
    clear Ttot_temp Ti_temp
    VFlow_crossing_i_saved = VFlow_crossing_i;
    VFlow_crossingE_i_saved = VFlow_crossingE_i;
    
    
    % This section looks at periods with gaps in breaths, and fills in with
    % smaller breaths: So nees to query whether should actually be dividing
    % into
    %[VFlow_crossing_i; VFlow_crossingE_i]
    
    for i=(length(VFlow_crossing_i)-2):-1:2
        Ttot_temp=time(VFlow_crossing_i(i+1))-time(VFlow_crossing_i(i));
        % Calculate a local median Ttot from the 5 previous, and 5
        % following periods, and check how it compares with the median
        % t-tot median for the whole period. If it is +/-35% of the median
        % value for the whole window, apply this to pad missing breaths:
        LocalRange=5;
        index_s=i-LocalRange; index_s(index_s<1)=1;
        index_e=i+LocalRange; index_e(index_e>length(VFlow_crossing_i)-1)=length(VFlow_crossing_i)-1;
        ttot_temp=diff(time(VFlow_crossing_i));
        LocalmedianTtot1=median(ttot_temp(index_s:index_e));
        medianTtot1;
        if (LocalmedianTtot1<0.65*medianTtot1)||(LocalmedianTtot1>1.35*medianTtot1)
            LocalmedianTtot1=medianTtot1;
        end
       
        % If there is a particularly long break (longer than a reasonable
        % apnea with no breath movement, treat as unsalvageable artefact,
        % and break:
        if Ttot_temp>20
            time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN;time=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; VFlow_filtered1=NaN; index=NaN;
            return;
        end
        
        if Ttot_temp>(2*LocalmedianTtot1)
            %add new crossings to make into 1 breath into X smaller breaths
            NnewBreaths = round(Ttot_temp/LocalmedianTtot1)-1;
            newTtots_i=round(Ttot_temp/(NnewBreaths+1)/dt);
            
            newi = VFlow_crossing_i(i):newTtots_i:(VFlow_crossing_i(i)+newTtots_i*NnewBreaths);
            newEi = newi+round(newTtots_i/2);
            
%             newEi = VFlow_crossingE_i(i):newTtots_i:(VFlow_crossingE_i(i)+newTtots_i*NnewBreaths) % This is problematic if newEi is not evenly distributed within the gaps - i.e. mid may be detected in breath.
%             VFlow_crossing_i = [VFlow_crossing_i(1:(i-1)) newi VFlow_crossing_i(i+2:end)];
%             VFlow_crossingE_i = [VFlow_crossingE_i(1:(i-1)) newEi VFlow_crossingE_i(i+2:end)];
            VFlow_crossing_i = [VFlow_crossing_i(1:(i-1)) newi VFlow_crossing_i(i+1:end)];
            VFlow_crossingE_i = [VFlow_crossingE_i(1:(i-1)) newEi VFlow_crossingE_i(i+1:end)];

            Ttot_temp(i)=time(VFlow_crossing_i(i+1))-time(VFlow_crossing_i(i));
        end
    end

end

%[VFlow_crossing_i; VFlow_crossingE_i]
% If less than 4 breaths detected in the period, set as artefact period. 
if length(VFlow_crossing_i)<4
    time1=NaN; VI1=NaN; BB_i_start=NaN; BB_i_mid=NaN; BB_i_end=NaN; BB_t=NaN; VTi=NaN; VTe=NaN; VI=NaN; VE=NaN; Ttot=NaN; Ti=NaN; Te=NaN; FLindex1=NaN; PeakInspFlow=NaN; PeakInspFlow_t=NaN; MidInspFlow=NaN; time=NaN; VFlow_filtered1=NaN; index=NaN;
    return;
end

%% Breath-breath data

clear Ttot BB_t BB_i Ti Te BB_i_mid BB_i_start BB_i_end VTi VTe;
Ti=zeros(1,length(VFlow_crossing_i)-1);
BB_i_mid=zeros(1,length(VFlow_crossing_i)-1);
Ttot=zeros(1,length(VFlow_crossing_i)-1);
BB_t=zeros(1,length(VFlow_crossing_i)-1);
BB_i_start=zeros(1,length(VFlow_crossing_i)-1);
BB_i_end=zeros(1,length(VFlow_crossing_i)-1);
VTi=zeros(1,length(VFlow_crossing_i)-1);
VTe=zeros(1,length(VFlow_crossing_i)-1);
PeakInspFlow=zeros(1,length(VFlow_crossing_i)-1);
PeakInspFlow_t=zeros(1,length(VFlow_crossing_i)-1);
MidInspFlow=zeros(1,length(VFlow_crossing_i)-1);
a_value=zeros(1,length(VFlow_crossing_i)-1);
FLindex1=zeros(1,length(VFlow_crossing_i)-1);
for i=1:length(VFlow_crossing_i)
    if i<length(VFlow_crossing_i)
        Ti(i)=time(VFlow_crossingE_i(i))-time(VFlow_crossing_i(i));
        BB_i_mid(i)=VFlow_crossingE_i(i);
        Ttot(i)=time(VFlow_crossing_i(i+1))-time(VFlow_crossing_i(i));
        BB_t(i)=time(VFlow_crossing_i(i));
        BB_i_start(i)=VFlow_crossing_i(i);
        BB_i_end(i)=VFlow_crossing_i(i+1)-1;
        VTi(i)=dt*sum(VFlow_filtered1_pos(VFlow_crossing_i(i):(VFlow_crossing_i(i+1)-1)));
        VTe(i)=-dt*sum(VFlow_filtered1_neg(VFlow_crossing_i(i):(VFlow_crossing_i(i+1)-1)));
        %prior to midinsponly:
        [PeakInspFlow(i),tempi]=max(VFlow_filtered1(VFlow_crossing_i(i):(VFlow_crossing_i(i)+round(0.5*Ti(i)/dt))));
        PeakInspFlow_t(i)=dt*(tempi-1);
        MidInspFlow(i)=VFlow_filtered1(VFlow_crossing_i(i)+round(0.5*Ti(i)/dt));
        a_value(i)=PeakInspFlow(i)./sin(pi*PeakInspFlow_t(i)./Ti(i));
        FLindex1(i)=(a_value(i)-MidInspFlow(i))/a_value(i);
    end
end
Te=Ttot-Ti;
VE = VTe./Ttot;
VI = VTi./Ttot;

FLindex1(FLindex1>1|FLindex1<0|isnan(FLindex1))=1;
FLindex2=1-MidInspFlow./PeakInspFlow; FLindex2(FLindex2>1|FLindex2<0|isnan(FLindex2))=1;

if ~minimum_figs;
figure(35);
Nplots=1;
ax35(1)=subplot(Nplots,1,1);
plot(time,VFlow_filtered1/max(VFlow_filtered1),VFlow_crossing_t,VFlow_crossing/max(VFlow_filtered1),'.r',VFlow_crossingE_t,VFlow_crossingE/max(VFlow_filtered1),'.k')
hold('on')
stairs(BB_t,FLindex1,'g');
stairs(BB_t,Ti./Ttot,'r');
stairs(BB_t,FLindex2,'k');
hold('off')
end

%% Make T=1 data
dT=1;
startT=ceil(BB_t(1));
endT=floor(BB_t(end));
time1 = startT:dT:endT;
% xmn=time(round((BB_i_start+BB_i_mid)/2))
% xmn(2:end)-xmn(1:end-1)
% round((BB_i_start+BB_i_mid)/2)
% BB_i_start
% BB_i_mid

VI1 = interp1(time(round((BB_i_start+BB_i_mid)/2)),VI,time1,'linear','extrap'); %VI is located at mid inspiration
VI1(VI1<0)=0;
% 
BB_i_start=BB_i_start+indexstart-1;
BB_i_mid=BB_i_mid+indexstart-1;
BB_i_end=BB_i_end+indexstart-1;

index=[indexstart indexend];
