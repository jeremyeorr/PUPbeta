%% Ventilatory control model validation
%
% Written by Philip I Terrill and Scott A Sands, 
% with the exception of the "pink noise function" provided by Hristo Zhivomirov, see below
% This MATLAB routine accompanies the manuscript entitled:
% "Quantifying the Ventilatory Control Contribution to Obstructive Sleep Apnea"

function [LGplusinfo_xx,LG1_xx,BiasLG1,AbsErr95]=ModelSimulation()
N=100;                        % choose number if simulations to run
LGcrit_xx=1.43*(rand(1,N)); % choose LG values to use
for xx=1:N
    %LGcrit_xx(xx)=1.43*(xx/N)+0.005*(rand(1,1)-0.5); % choose LG values to use
    [LGplusinfo_xx(xx,:),LG1_xx(xx),gamma1(xx),LGSS_xx(xx),LGn_xx(xx),LGcrit_xx(xx)]=RunTheModel(LGcrit_xx(xx),xx);  % Provides equally spaced samples up to a LG1 of 2)
    figure(2);
    plot(LG1_xx,LGplusinfo_xx(:,6),'MarkerSize',14,'Marker','.','LineStyle','none','Color',[0.2 0.1 0.65]);
    set(gcf,'color',[1 1 1])
    set(gca,'box','off','fontname','arial narrow')
    xlabel('LG1 true'); ylabel('LG1 estimated');
    BiasLG1=mean(LGplusinfo_xx(:,6)'-LG1_xx);
    AbsErr95=prctile(abs(LGplusinfo_xx(:,6)'-LG1_xx),95)
end
end

%% Model Simulation

function [LGplusinfo,LG1,gamma1,G_,LGn,LGcrit]=RunTheModel(LGcrit,xx)
%This function simulates OSA for a single value of loop gain, LGcrit:
WindowSize=7; % 7 minute analysis windows
plot_vent=1;
prop=25;  %determines approximate proportion of time in events, %
amplitude=0.01; %relative noise amplitude (SD/mean), Fr
%**********************************************************************
% Define The model parameters for simulated ventilation, and set
% initial values of key simulation variables
%**********************************************************************
dt=0.01; %Sampling interval
Ttot=3.5;    %Breath interval
a=-Ttot*2; %Start time
b=WindowSize*60+dt; %End time
time = a:dt:b;  %Time array
N=length(time); %No. of samples

VL = 3.5;                               %lung volume, L
PI_CO2_ = 0;                            %inspired carbon dioxide, mmHg
KCO2 = 0.0048;                          %blood capacitance for CO2
MRCO2 = 0.2/60;                         %metabolic production of CO2, L/s
Qdotp = 3/60;                           %Cardiac output / pulmonary blood flow, L/s
Pa_CO2_ = 40;                           %Resting/eupneic alveolar pCO2, mmHg
delay = 12;                             %circulatory delay
tau = 0;                                %time constant of ventilatory response to CO2 (default=0, instantaneous)
VE_ = 863*MRCO2/(Pa_CO2_-PI_CO2_);   %find resting ventilation

%The critical chemoreflex gain that just results in instability (Gcrit) depends on the delay. List is relevant for default parameters only
switch delay
    case 20,	Gcrit = 0.011881178787224/1.000000000000007; %(tau=0, delay=20)
    case 25,	Gcrit = 0.010619083081828/0.999999999999973; %(tau=0, delay=25)
    case 24,    Gcrit = 0.010826525587416/1.000000000000027; %(tau=0, delay=8)
    case 16,    Gcrit = 0.013501440565352/0.999999999999997; %(tau=0, delay=8)
    case 12,	Gcrit = 0.016259528066649/1.000000000000012; %(tau=0, delay=12)
    case 8,     Gcrit = 0.021873092173401/0.999999999999990; %(tau=0, delay=8)
    case 6,     Gcrit = 0.027540222444834/1.000000000000015; %(tau=0, delay=6)
    case 4,     Gcrit = 0.038932124036743/1.000000000000002; %(tau=0, delay=4)
end
G = LGcrit*Gcrit;

Ip = Pa_CO2_ - VE_/G;                    % Set threhold CO2 for apnea based on the gain
tau_washout = VL/(VE_+Qdotp*KCO2*863);   % Plant time constant
Gplant=Pa_CO2_/VL;                       % Plant gain term constant
G_=G*Gplant*tau_washout;

% Subtract Ttot from delay when using sample&hold simulations since breath time adds additional delay with respect to instability.
delayX = delay-Ttot/2; %discrete nature of ventilation adds an equivalent delay of between 0-Ttot
delay_i = round(delayX/dt);

% Initialize continuous variables
Pa_CO2 = zeros(1,N); Pa_CO2(1)=40;
Pa_CO2_p = zeros(1,N);
VE = zeros(1,N);
VE_chem = zeros(1,N);
VE_output = zeros(1,N);
%**********************************************************************
%**********************************************************************

%**********************************************************************
% Using the theoretical ventilatory model, and the defined parameters
% calculate the actual loop gain variables.
%**********************************************************************
F=[0.002:0.00001:1.09999 1.1.^(1:100)];     %frequency, Hz
phase=pi-atan(2*pi*tau*F)-atan(2*pi*tau_washout*F)-2*pi*delay*F;
F180=interp1(phase,F,0); % natural frequency
T180=1/F180;       %natural period

%Calculate true loop gain:
LGn=abs(-G*Gplant*exp(-1i*2*pi*F180*delay)/(1/tau_washout+1i*2*pi*F180)/(1+1i*2*pi*F180*tau));
%LGn=loop gain at the natural oscillatory frequency (LGn>1 is unstable)
LG1=abs(-G*Gplant*exp(-1i*2*pi*(1/60)*delay)/(1/tau_washout+1i*2*pi*(1/60))/(1+1i*2*pi*(1/60)*tau))
%LG1=loop gain at 1 cycle/min i.e. heart failure
LG2=abs(-G*Gplant*exp(-1i*2*pi*(1/30)*delay)/(1/tau_washout+1i*2*pi*(1/30))/(1+1i*2*pi*(1/30)*tau));
%LG2=loop gain at 2 cycle/min i.e. idiopathic central sleep apnea
LGf=(-G*Gplant*exp(-1i*2*pi*F*delay)./(1/tau_washout+1i*2*pi*F)./(1+1i*2*pi*F*tau));
%Full loop gain vs frequency function

%**********************************************************************
%**********************************************************************


%**********************************************************************
% In order to simulate ventilation in OSA patients disturbance is added
% to the theoretical ventilatory control model:
% - Random white noise added to ventilatory drive
% - Obstructive respiratory events simulated with a gradual reduction
% of flow
% - An additional ventilatory response to arousal which occurs with a
% defined probability at the end of a respiratory event
% - An additional ventilatory response to arousal which occurs with a
% defined probability at any breath
%**********************************************************************

%Ventilatory response to arousal (0.4 = 40% extra ventilation above the resting level)
gamma1=0.4; %Set this.

%Probability of arousals at end of event and spontaneously
Ar_chance1=0.8; Ar_chance2=0.005;

%Make a random ventilatory (white noise) disturbance
M=ceil(N/(Ttot/dt));

time_ds = a:Ttot:b;
VErand = VE_*(randn(1,M));
VErand=VErand/(std(VErand)/VE_); %correction for precise power.
VE_disturbance=interp1(time_ds,amplitude*VErand,time,'nearest','extrap'); %could use nearest

%Create randomly timed respiratory events: pink noise works quite well
temp = detrend(pinknoise(M));
%
thres = prctile(temp,prop);
RDist = zeros(1,M);
Ar = zeros(1,M);
VDist = zeros(1,M);
for i=3:length(temp)    %Events occur after a certain time if below a threshold for 3 breaths
    if time_ds(i)>0&&(temp(i)<thres||temp(i-1)<thres||temp(i-2)<thres)
        RDist(i)=1;
    end
end

%Apply arousals randomly:
for i=2:M
    if (RDist(i)==0)&&(RDist(i-1)==1)&&(rand(1)<Ar_chance1) %resp arousal
        Ar(i)=1;
    end
    if (RDist(i)==0)&&(rand(1)<Ar_chance2) %spontaneous arousal
        Ar(i)=1;
    end
end

%Make arousals last for two breaths.
tempAr=Ar;
if 1
    for i=1:length(Ar)-1
        if Ar(i)==1
            tempAr(i+1)=1;
        end
    end
end
Ar=tempAr;

%Remove resistance/events at arousal onset for 2 breaths
RDist([Ar 0]>0)=0;

%gradual reduction in flow with obstruction...
shape=[0.25 0.5 0.75];
for i=2:length(RDist)
    if RDist(i)==1&&RDist(i-1)==0
        RDist(i)=shape(1);
    end
    if RDist(i)==1&&RDist(i-1)==shape(1)
        RDist(i)=shape(2);
    end
    if RDist(i)==1&&RDist(i-1)==shape(2)
        RDist(i)=shape(3);
    end
end
Obstruction=interp1(time_ds,RDist,time,'nearest','extrap');

for i=2:M
    if Ar(i)==1
        VDist(i)=VE_*gamma1;
    end
end
Ar_us=interp1(time_ds,Ar,time,'nearest','extrap');
VE_disturbance=VE_disturbance+interp1(time_ds,VDist,time,'nearest','extrap');
%**********************************************************************
%**********************************************************************

%**********************************************************************
% Observed ventilation does not have a constant Ttot - rather, Ttot
% varies from breath to breath. In order to simulate this in model
% derived ventilation, a mean Ttot is defined, with a randomised
% variability.
%**********************************************************************
time_ds = a:Ttot:b;
time_ds = sort((time_ds + 0.05*Ttot*randn(1,length(time_ds))));
time_ds(1)=a;
time_ds(end)=b;
BreathTime = 0*time;
for i=2:N
    for x=1:length(time_ds)
        if time(i)>=time_ds(x)&&time(i-1)<time_ds(x)
            BreathTime(i)=1;
        end
    end
end
%**********************************************************************
%**********************************************************************

%**********************************************************************
% Generate simulated ventilation. Model parameters, programmed
% ventilatory disturbance, and the defined variable Ttot series.
%**********************************************************************
for i = 1:N,
    
    % Handle initial conditions for delayed ventilation in inital
    % segment of the model
    if (i-delay_i)>=1
        if tau>0
            Pa_CO2_p(i)  = Pa_CO2_p(i-1) + dt/tau*(Pa_CO2(i-delay_i) - Pa_CO2_p(i-1));
        else
            Pa_CO2_p(i)  = Pa_CO2(i-delay_i);
        end
    else
        Pa_CO2_p(i)  = Pa_CO2(i);
    end
    
    VE_chem(i) = G*(Pa_CO2_p(i)-Ip); %Linear chemoreflex response
    VE_output(i) = VE_chem(i) + VE_disturbance(i); % Add random noise (as described above) to chemical drive
    
    %If resistance increases ventilatory output is redued
    if Obstruction(i)>0
        VE_output(i) = VE_output(i)*(1-0.5*Obstruction(i));
    end
    
    % Ventilation is zero if ventilatory drive is negative
    if VE_output(i)<0
        VE(i)=0;
    else
        VE(i)=VE_output(i);
    end
    
    if i==1
        currentV=VE(i);
    end
    
    % Instantaneous level of VDr is kept for the whole breath
    if BreathTime(i)
        currentV=VE(i);
    else
        VE(i)=currentV;
    end
    
    % Iterate through such that final element of Pa_CO2 is the Nth
    % element
    if i<N
        Pa_CO2(i+1)=Pa_CO2(i)+dt*(-Gplant*(VE(i)-VE_)-(Pa_CO2(i)-Pa_CO2_)/tau_washout);
    end
    
end
%**********************************************************************
%**********************************************************************

%**********************************************************************
% Use the simulated ventilation and simulated respiratory events as
% inputs to validate the estimation of loop gain parameters using the
% PSG loop gain estimation method
%**********************************************************************
VI=VE(BreathTime==1);
TimeB=time(BreathTime==1);
Vchem_true=VE_chem(BreathTime==1);
Ar=Ar_us(BreathTime==1);

CentralApneas=0*VI; CentralApneas(Vchem_true<0)=1;

ObstructionB=Obstruction(BreathTime==1);
ObstructedBreaths=1+0*VI; ObstructedBreaths(ObstructionB>0)=0;
ObstructedBreaths(CentralApneas==1)=0;

UseTrueEupnea=1; %In practise, we do not know true eupnea. Thus we use the mean.
if UseTrueEupnea
    Veupnea=VE_;
else
    Veupnea=mean(VI');
end

Ttot_n=[Ttot diff(time(BreathTime==1))];
Ttot_previous=[diff(time(BreathTime==1)) Ttot];

%     [Breath-by-Breath minute ventilation      Obstructed Breaths     Arousal Breaths     The time of the breath  Central apnoea breaths   Ttot series for current breath   Ttot series for previous breath  ]
xdata=[VI'-Veupnea                           ObstructedBreaths'     Ar'                 TimeB'                  CentralApneas'           Ttot_n'                        Ttot_previous'                   ];
VarousalOn=1; % VarousalOn=1 to fit a ventilatory response to arousal parameter to observed ventilation
polyfitorder=3;
[Error,Vchem_est,Varousal_est,LGplusinfo,BestParameters,BestSSres,FitQuality,i_1,i_end] = FindBestModelParameters(xdata,VarousalOn,Veupnea,polyfitorder); % This function fits a steady-state gain, time constant, delay and VRA model to observed ventilation to estimate loop gain
%**********************************************************************
%**********************************************************************


%**********************************************************************
% Plot the simulated ventilation, and the model fit to the "observed"
% ventilation.
%**********************************************************************
if plot_vent==1
    plot_multiple=1; %set to '1' to display multiple figures
    gplot=1;
    figure(gplot);
    
    if plot_multiple
        Nsubplots=8;
        cols=2;
        nsubplot=1+mod(xx-1,Nsubplots);
        ax1(nsubplot)=subplot(Nsubplots/cols,cols,nsubplot);
        fontsize_=9;
    else
        ax1(1)=subplot(1,1,1);
        fontsize_=12;
    end
    
    PlotModelResults(LGplusinfo,xdata,Vchem_est,Varousal_est,Veupnea,Vchem_true,LG1,i_1,fontsize_)
end

end

%% Plot Results Function
function PlotModelResults(LGplusinfo,xdata,Vchem_est,Varousal_est,Veupnea,Vchem_true,LG1,i_1,fontsize_)
E1=xdata(:,2);
CentralApneas=xdata(:,5);
Ar=xdata(:,3);
BB_t=xdata(:,4);
%BB_t is the start time of each breath
%BB_i_start is the start index (i.e. row) of each breath

set(gcf,'Color',[1 1 1]);
dt_plot=0.25;
startT=ceil(BB_t(1));
endT=ceil(BB_t(length(BB_t)));
time_dt=startT:dt_plot:endT;
E1_rs=0*time_dt;
for i=1:length(time_dt);
    E1_rs(i) = 1-E1(find(BB_t<=time_dt(i),1,'last'));
end
C1_rs=0*time_dt;
for i=1:length(time_dt);
    C1_rs(i) = CentralApneas(find(BB_t<=time_dt(i),1,'last'));
end
E1_rs=E1_rs&(~C1_rs);
Ar_rs=0*time_dt;
for i=1:length(time_dt);
    Ar_rs(i) = 1-Ar(find(BB_t<=time_dt(i),1,'last'));
end

tempsize=0.25; tempoffset=0.375;
tempy=tempsize*(E1_rs)+tempoffset;

hold('off')
%plot obstructive event breaths over the top of flow in blue
tempsize=0.5; tempoffset=0.75;
tempy=tempsize*(E1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
hold('on');
%plot central apnea event breaths over the top of flow in red
tempsize=0.5; tempoffset=0.75;
tempy=tempsize*(C1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);

tempoffset=1.1*max((xdata(1:end,1)+Veupnea)/Veupnea); tempsize=tempoffset/15;
tempy=tempsize*(1-Ar_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[0.1 1 0.1],'EdgeColor','none','FaceAlpha',0.5);

if ~isempty(LG1)
    ylabel(['True LG=' num2str(LG1,2) ', Est. LG=' num2str(LGplusinfo(6),2)],'fontsize',fontsize_,'fontname','arial narrow')
else
    ylabel(['Est. LG=' num2str(LGplusinfo(6),2)],'fontsize',fontsize_,'fontname','arial narrow')
end
hold('on');
stairs(BB_t(1:end),(xdata(1:end,1)+Veupnea)/Veupnea,'color',[0.4 0.4 0.7]);
if ~isempty(Vchem_true)
    plot(BB_t(i_1:end),(Vchem_true(i_1:end)')/Veupnea,'--','color',[0 0 0]);
end
stairs(BB_t(i_1:end)+0.3,(Vchem_est(i_1:end)'+Varousal_est(i_1:end)'+Veupnea)/Veupnea,'color',[0.1 0.8 0.1]);
plot(BB_t(i_1:end),(Vchem_est(i_1:end)'+Veupnea)/Veupnea,'color',[0 0 0]);

set(gca,'FontSize',fontsize_,'FontName','Arial Narrow','tickdir','out');
box('off');
xlim([min(BB_t) max(BB_t)+5]);
hold('off');

Ylim1=get(gca,'Ylim')
if Ylim1(1)>0
    Ylim1(1)=0;
end
set(gca,'Ylim',Ylim1)
end

%% Pink noise function
% see http://www.mathworks.com/matlabcentral/fileexchange/42919-pink-red-blue-and-violet-noise-generation-with-matlab-implementation/content/pinknoise.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pink Noise Generation with MATLAB Implementation   %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov       07/30/13  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = pinknoise(N)
% function: y = pinknoise(N)
% N - number of samples to be returned in row vector
% y - row vector of pink (flicker) noise samples
% The function generates a sequence of pink (flicker) noise samples.
% Pink noise has equal energy in all octaves (or similar log bundles) of frequency.
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB per octave.
% define the length of the vector
% ensure that the M is even
if rem(N,2)
    M = N+1;
else
    M = N;
end
% generate white noise with sigma = 1, mu = 0
x = randn(1, M);
% FFT
X = fft(x);
% prepare a vector for 1/f multiplication
NumUniquePts = M/2 + 1;
n = 1:NumUniquePts;
n = sqrt(n);
% multiplicate the left half of the spectrum so the power spectral density
% is inversely proportional to the frequency by factor 1/f, i.e. the
% amplitudes are inversely proportional to 1/sqrt(f)
X(1:NumUniquePts) = X(1:NumUniquePts)./n;
% prepare a right half of the spectrum - a copy of the left one,
% except the DC component and Nyquist frequency - they are unique
X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));
% IFFT
y = ifft(X);
% prepare output vector y
y = real(y(1, 1:N));
% normalise
y = y./max(abs(y));
end