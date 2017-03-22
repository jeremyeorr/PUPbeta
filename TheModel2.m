function [Vchem] = TheModel2(Parameters,Data,i_1,delayedVE,polyfitorder)
% Written by Philip I Terrill and Scott A Sands
% This MATLAB function accompanies the manuscript entitled:
% "Quantifying the Ventilatory Control Contribution to Obstructive Sleep Apnea"
%
% This function is used to identify the error between the
% predicted ventilatory drive for a given set of parameters (LG0, Tau, gamma
% delay, Error[0]), and the observed ventilation when the airway is open. 
% Mean square error and the r_squared value is provided to quantify the quality of the model fit.
% This function is potentially iterated thousands of times by the fmincon
% function call within the function "Ventilation_Model_Fit.m". As such all possible processing is
% conducted in other parts of the program.
%
% >>>>>>>>>>>>>>>>>>>>>>>>Input parameters<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% - Parameters: This is an array 1*4 array of the parameters to be fit: [LG0         tau         Error[0]          gamma]
% - Data=[VI'-Veupnea E1' AR1(1:length(E1))' TimeB' Ecentralapnea' [Ttot_ser.'; 0]];
% ***xdata: This matrix contains information about the patients observed
% ventilation
%             - Column 1: Breath-by-Breath Minute Ventilation
%             - Column 2: Obstructive respiratory event series. Elements=1 for an un-obstructed breath; Elements=0 for an obstructed breath.
%             - Column 3: Arousal Scoring series: defines whether there is an arousal present in each breath (arousal present, AR[n]=1; No arousal present, AR[n]=0);
%             - Column 4: The time that each breath in the series occurs. 
%             - Column 5: Scored Central Apnea array series: defines whether each breath is scored as a central apnea. elements=1 for breaths in a scored central apnea; Elements=0 for breaths not classified as a central apnea.
%             - Column 6: T_tot series (unused here)
%             - Column 7: T_tot series shifted by 1 breath (previous Ttot)
% ***i_1: This is the index of the first breath in the series (contained in
% xdata) that the model is fit to.
% ***delayedVE: This is the series of delayed ventilation data (generated from the observed ventilation data using a linear interpolation method between breaths)
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
% >>>>>>>>>>>>>>>>>>>>>>>>Output parameters<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% - F: 1*1 Array. Mean Square Error (MSE) of the model fit to observed
% ventilation data
% - Error: 1*N Array. This is the breath-by-breath different between the
% model predicted ventilation and the observed ventilation
% - Rsq: The R-squared value of the linear correlation between the model
% estimated ventilation and the observed ventilation, including the
% polynomial drift correction
% - Vchem_est: The ventilatory drive due to chemical drive
% - VAr_est: The ventilatory drive due to arousal
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

VE=Data(:,1);               %The observed breath-by-breath minute ventilation
Ttot_previous=Data(:,7);    %The previous breath duration (Ttot)
LG0=Parameters(1);          %Estimated model parameters
tau1=Parameters(2);
Err_0=Parameters(3);
gamma=Parameters(4); 

%**************************************************************************
% Ventilatory response to arousal is modelled using the additional
% parameter gamma, that defines additional ventilatory drive in breaths
% that are temporally associated with an arousal. i.e. For each breath n we
% determine whether there is a scored arousal present, Ar[n]=1, or not,
% Ar[n]=0, such that Varousal is simply given by Equation S4:
% Varousal[n]=Gamma*Ar[n].
%**************************************************************************
AR=Data(:,3);               %This array defines whether there is an arousal present in each breath (arousal present, AR[n]=1; No arousal present, AR[n]=0);
Varousal=gamma*AR;          %Estimated Varousal trace
%**************************************************************************
%**************************************************************************

%**************************************************************************
% Calculate the compound variables alpha and beta. As per formula E1:
% Vchem[n]=alpha*Vchem[n-1]+beta*VE'[n], where alpha and beta are calculated using a
% forward and backward euler approximation. 
% Note, alpha and beta vary with each breath according to changes in T_tot[n-1].
% As such, these variables are vectors of the same length as observed
% ventilation. 
%**************************************************************************
alpha1=(tau1./Ttot_previous)./(1+tau1./Ttot_previous); % Forward Euler approximation of alph
alpha2=1-Ttot_previous/tau1; % Backward Euler approximation of alpha
beta1=-LG0./(1+tau1./Ttot_previous); % Forward Euler approximation of beta
beta2=-LG0./(tau1./Ttot_previous); % Backward Euler approximation of beta
alpha=(alpha1+alpha2)/2; % alpha
beta=(beta1+beta2)/2; % beta
%**************************************************************************
%**************************************************************************

%**************************************************************************
% Iterate through the ventilation series to calculate the chemical
% ventilatory drive (V_chem) at each breath using eqution E1:
% Vchem[n]=alpha*Vchem[n-1]+beta*delayedVE[n]
%**************************************************************************
Vchem(i_1-1)=VE(i_1-1)-Varousal(i_1-1)-Err_0; % Chemical drive must be seeded with an initial value of chemical drive at time=0 (breath n=1), using the model parameter Vchem[1]

for k=i_1:length(Data)
    Vchem(k)=alpha(k)*Vchem(k-1)+delayedVE(k)*beta(k);
end
%**************************************************************************
%**************************************************************************

%**************************************************************************
% Overall ventilatory drive is modelled as the combined effect of chemical
% drive and ventilatory response to arousal (Equation 1):
% Vdrive=Vchem+Varousal	(Equation 1)
%**************************************************************************
Vdrive=Vchem.'+Varousal;
%**************************************************************************
%**************************************************************************

%**************************************************************************
% It is now possible to calculate the error  between the predicted
% ventilation and the and the observed ventilation (as per Equation E6):
% Error[n]=Vdrive[n]-VE[n]
%**************************************************************************
Error1=VE-Vdrive;
%**************************************************************************
%**************************************************************************


    %**************************************************************************
    % We now aim to determine how well the model predicted ventilation fits the
    % observed ventilation. The model fit is quantified by the
    % mean-square-error. However, we consider three physiological factors in
    % this calculation:
    %
    % 1. During obstruction, ventilatory drive is not
    % expressed in observed ventilation. As such we penalize only
    % non-obstructed breaths.
    %
    % 2. During a central apnea, however, observed
    % ventilation is ~zero, but the Vdrive is expected to fall below
    % zero (e.g. CO2 below the apneic threshold); during these times we do not
    % expect VE to reflect Vdrive. Thus the error is instead taken as 0 if
    % estimated Vdrive is below 0 (thereby avoiding an unreasonable penalty)
    %
    % 3. To deal with non-stationarity (e.g. drifting signal amplitude) and
    % non-Gaussian (non-white) error, we subtracted a third-order
    % polynomial P[n] from the model error (Equation S6), thereby removing the
    % effect of long-term error or drift: Error*[n]=Error[n]xP[n].
    %**************************************************************************
    Penalize=Data(:,2); % The Penalize array, also known as W[n], defines whether each breath is obstructed. Penalize[n]=1 for an un-obstructed breath; Penalize[n]=0 for an obstructed breath.
    CentralApneas=Data(:,5); % The CentralApneas array defines whether each breath is scored as a central apnea. CentralApneas[n]=1 for a central apnea; CentralApneas[n]=0 for a breath not classified as a central apnea.
    Penalize(CentralApneas==1&Error1>0)=0; % Remove penalties for central apneas with positive noise (Drive<0)
    Penalize(Penalize==0&CentralApneas==0&Error1>0)=1; %positive noise greater than tolerance is not expected during OBstr
    Penalized_Error1=Error1.*Penalize; % The array Penalised_Error defines the error for each penalised (unobstructed) breath.
    %**************************************************************************
    % A third order polynomial if fit to the derived breath-by-breath error,
    % and this polynomial trend is removed from the error signal.
    %**************************************************************************
    tempW=Error1(i_1:end);
    tempE1=Penalize(i_1:end);
    tempy=Data(:,4);
    tempy2=Data((i_1:end),4);
    polyfitresultsErr=polyfit(tempy2(tempE1>0),tempW(tempE1>0),polyfitorder);
    ErrMovingFit=0*tempy;
    ErrMovingFit(i_1:end)=polyval(polyfitresultsErr,tempy(i_1:end));
    
    %Penalized_Error=0*Penalized_Error1;
    if 1 %Subtract the polynomial drift, but leaving the 'mean' error untouched. Solutions whose mean is non-zero will not be favored.
    Penalized_Error(Penalize>0)=Penalized_Error1(Penalize>0)-ErrMovingFit(Penalize>0)+mean(ErrMovingFit(i_1:end));
    Error=Error1-ErrMovingFit+mean(ErrMovingFit(i_1:end));
    Vchem=Vchem.'+ErrMovingFit-mean(ErrMovingFit(i_1:end));
    else %Option to not use the "mean subtraction" procedure.
    Penalized_Error(Penalize>0)=Penalized_Error1(Penalize>0)-ErrMovingFit(Penalize>0);
    Error=Error1-ErrMovingFit;
    Vchem=Vchem.'+ErrMovingFit;
    end
    
    Error(1:i_1-1)=NaN; % Set values before breath n=1 (i_1) to unknown or "NaN", since these are not fitted values
    Vchem(1:i_1-1)=NaN; % Again, set values before i_1 to NaN, since these are not fitted values
    %**************************************************************************
    % Based on the filtered error signal, calculate the metrics to describe how
    % well the model fits the observed ventilation data
    %**************************************************************************
    SSres=(sum(Penalized_Error(i_1:end).^2))/(sum(Penalize(i_1:end))); % Final Mean Square Error for unobstructed breaths.
        
    SStot_temp=VE(Penalize>0)-mean(VE(Penalize>0));
    SStot=sum(SStot_temp(i_1:end).^2)/(sum(Penalize(i_1:end)));
    Rsq=1-SSres/SStot;%(sum(Penalized_Error(i_1:end).^2))/(sum(SStot(i_1:end).^2)); % The R_squared value for how well the predicted ventilation correlates with the observed ventilation data.
    
    
    
    %**************************************************************************
    %**************************************************************************

end


