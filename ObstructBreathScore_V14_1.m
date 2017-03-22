function [E1,Ecentralapnea,Ecentralhypop] = ObstructBreathScore_V14_1( Veupnea, E, VI, VEThorax, VEAbdo, BB_i_start, removecentralhypops, findcentralhypopneasandapneas,seteventsbasedonvlessthaneupnea2,eventslessthanX, plot_on)


    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Update events based on eupnea, and detection of "detected" central
    % apneas... i.e. any breath greater than eupnoea has events cleared,
    % any breath at less than 70% of eupnoea is set as an event. Also, any
    % breath detected as a central event (depending on criteria) has event
    % cleared for purpose of model fitting. 
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    E1=E;
    EventBreaths=VI;
    EventBreaths(E==1)=[];
    ClearBreaths=VI;
    ClearBreaths(E==0)=[];

   % Clear manual events, and apply a hard threshold to identify events.  
   if seteventsbasedonvlessthaneupnea2==1
        E1(1:end)=1;
        for k=1:length(BB_i_start)
            if (VI(k)<Veupnea*eventslessthanX)
                E1(k)=0;
            end
        end
   % Maintain manual events, but edit to score any below hard threshold        
   elseif seteventsbasedonvlessthaneupnea2==2
        for k=1:length(BB_i_start)
            if (VI(k)<Veupnea*eventslessthanX)
                E1(k)=0;
            elseif (VI(k)>Veupnea)
                E1(k)=1;
            end
        end
   % Maintain manual events, and use histogram to define bounds for re-classifiying breaths.      
   elseif seteventsbasedonvlessthaneupnea2==3
        for k=1:length(BB_i_start)
            if (VI(k)<max([Veupnea*eventslessthanX median(EventBreaths)]))
                E1(k)=0;
            elseif (VI(k)>median(ClearBreaths))
                E1(k)=1;
            end
        end
    
   elseif seteventsbasedonvlessthaneupnea2==4
        for k=1:length(BB_i_start)
            if (VI(k)<median(EventBreaths))
                E1(k)=0;
            elseif (VI(k)>median(ClearBreaths))
                E1(k)=1;
            end
        end
        
   elseif seteventsbasedonvlessthaneupnea2==5
       if ((1-(sum(E)/length(E)))>0.15)&&((1-(sum(E)/length(E)))<0.85)&&(median(EventBreaths)<median(ClearBreaths))
           for k=1:length(BB_i_start)
                if (VI(k)<median(EventBreaths))
                    E1(k)=0;
                elseif (VI(k)>median(ClearBreaths))
                    E1(k)=1;
                end
           end 
       else
           for k=1:length(BB_i_start)
               if (VI(k)<Veupnea*eventslessthanX)
                   E1(k)=0;
               elseif (VI(k)>Veupnea)
                   E1(k)=1;
               end
           end
       end
        
   % Only edit events at the edge of existing manually defined events:     
   elseif seteventsbasedonvlessthaneupnea2==6
        if ((1-(sum(E)/length(E)))>0.15)&&((1-(sum(E)/length(E)))<0.85)&&(median(EventBreaths)<median(ClearBreaths))
            ThresLow=median(EventBreaths);
            ThresHigh=median(ClearBreaths);
        else
            ThresLow=Veupnea*eventslessthanX;
            ThresHigh=Veupnea;
        end
       
        E1=E;
        E_Trans=E;
        E_Trans(1:end)=0;
        EventMod=1;
        while EventMod==1
            % Identify the breaths before and after events, check
            % threshold, and mark as events
            E_temp=E1;
            for k=1:length(BB_i_start)-1
                if (E1(k)==1)&&(E1(k+1)==0) % k is breath before event starts; k+1 is the first breath of event
                    E_Trans(k)=E_Trans(k)+1;
                    E_Trans(k+1)=E_Trans(k+1)+-1;
                elseif (E1(k)==0)&&(E1(k+1)==1) % k is last breath of event; k+1 is first breath after event
                    E_Trans(k)=E_Trans(k+1)+-1;
                    E_Trans(k+1)=E_Trans(k+1)+1;
                end    
            end   
            E1((VI<ThresLow)&(E_Trans>0))=0; % for breaths identified as the breath before or after scored event (even if it is a single non event breath) and ventilation below threshold, score as an event breath

            % Identify first and last breaths of event, check
            % threshold, and clear events as appropriate. 
            for k=1:length(BB_i_start)-1
                if (E1(k)==1)&&(E1(k+1)==0) % k is breath before event starts; k+1 is the first breath of event
                    E_Trans(k)=1;
                    E_Trans(k+1)=-1;
                elseif (E1(k)==0)&&(E1(k+1)==1) % k is last breath of event; k+1 is first breath after event
                    E_Trans(k)=-2;
                    E_Trans(k+1)=2;
                end    
            end   
            E1((VI>ThresHigh)&(E_Trans==-1))=1; % For first, or last breath of event, unless they are the same breath (when E_trans=-2, not -1), and ventilation greater than threshold, clear event.

            if sum(abs(E_temp-E1))==0 %i.e. no changes to events
                EventMod=0;
            else
                EventMod=1;
            end
        end
    % Only edit events at the edge of existing manually defined events, and
    % only using the >1 to clear, and <0.7 to mark. for between 0.7 and 1,
    % events left as-is
    elseif seteventsbasedonvlessthaneupnea2==7
        E1=E;
        E_Trans=E;
        E_Trans(1:end)=0;
        EventMod=1;
        while EventMod==1
            % Identify the breaths before and after events, check
            % threshold, and mark as events
            E_temp=E1;
            E_Trans(1:end)=0;
            for k=1:length(BB_i_start)-1
                if (E1(k)==1)&&(E1(k+1)==0) % k is breath before event starts; k+1 is the first breath of event
                    E_Trans(k)=1;
                    E_Trans(k+1)=-1;
                elseif (E1(k)==0)&&(E1(k+1)==1) % k is last breath of event; k+1 is first breath after event
                    E_Trans(k)=-1;
                    E_Trans(k+1)=1;
                end    
            end   
            E1((VI<Veupnea*eventslessthanX)&(E_Trans>0))=0; % for breaths identified as the breath before or after scored event (even if it is a single non event breath) and ventilation below threshold, score as an event breath

            % Identify first and last breaths of event, check
            % threshold, and clear events as appropriate. 
            E_Trans(1:end)=0;
            for k=1:length(BB_i_start)-1
                if (E1(k)==1)&&(E1(k+1)==0) % k is breath before event starts; k+1 is the first breath of event
                    E_Trans(k)=E_Trans(k);
                    E_Trans(k+1)=E_Trans(k)-1;
                elseif (E1(k)==0)&&(E1(k+1)==1) % k is last breath of event; k+1 is first breath after event
                    E_Trans(k)=E_Trans(k)-1;
                    E_Trans(k+1)=E_Trans(k)+1;
                end    
            end   
            E1((VI>Veupnea)&(E_Trans==-1))=1; % For first, or last breath of event, unless they are the same breath (when E_trans=-2, not -1), and ventilation greater than threshold, clear event.

            if sum(abs(E_temp-E1))==0 %i.e. no changes to events
                EventMod=0;
            else
                EventMod=1;
            end
        end
    end
   
    EventBreaths1=VI;
    EventBreaths1(E1==1)=[];
    ClearBreaths1=VI;
    ClearBreaths1(E1==0)=[];
    N_EB = hist(EventBreaths,0.05:0.10:2);
    N_CB = hist(ClearBreaths,0.05:0.10:2);
    N_EB1 = hist(EventBreaths1,0.05:0.10:2);
    N_CB1 = hist(ClearBreaths1,0.05:0.10:2);
    if plot_on==1
        figure
        subplot(1,2,1)
        hold on
        plot(0.05:0.10:2,N_EB,'r')
        plot(0.05:0.10:2,N_CB,'b')
        subplot(1,2,2)
        hold on
        plot(0.05:0.10:2,N_EB1,'r')
        plot(0.05:0.10:2,N_CB1,'b')
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Identify any central events based on the derived VE data, and
    % the abdo and thorax RIP data:
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Ecentralhypop=0*E;
    Ecentralapnea=0*E;
    if findcentralhypopneasandapneas
        TolCE=0.0;
        ThresCAVE=0.2;
        ThresCARIP=0.2;
        for k=1:length(BB_i_start)
            if VI(k)<1&&(E(k)||(VI(k)<eventslessthanX)) %If ventilatyion is lees than mean and there is a scored apnea/hypopnea or ventilation is <0.7:
                if (VI(k)<ThresCAVE)&&(VEThorax(k)<ThresCARIP)&&(VEAbdo(k)<ThresCARIP)
                    Ecentralapnea(k)=1;
                elseif (VI(k)>VEThorax(k)-TolCE)&&(VI(k)>VEAbdo(k)-TolCE) %If ventilation is less than 20% of eupnea, and thorax and abdomen are <30% of eupneic level, assume no drive, and therefore probably central event rather than a central event
                    Ecentralhypop(k)=1;
                end
            end

        end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Identify any single breath central's, and any single non-centrals
    % (within a central cluster) and remove. 
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    for k=2:(length(VI)-1)
        %if in a sequence of obstructed breaths... (three breaths of
        %event... k is middle breath)
        if ((E(k)==0)&&(E(k-1)==0)&&(E(k+1)==0))
            if Ecentralhypop(k)==1&&Ecentralhypop(k-1)==0&&Ecentralhypop(k+1)==0
                Ecentralhypop(k)=0;
            elseif Ecentralhypop(k)==0&&Ecentralhypop(k-1)==1&&Ecentralhypop(k+1)==1
                Ecentralhypop(k)=1;
            end
            if Ecentralapnea(k)==1&&Ecentralapnea(k-1)==0&&Ecentralapnea(k+1)==0
                Ecentralapnea(k)=0;
            elseif Ecentralapnea(k)==0&&Ecentralapnea(k-1)==1&&Ecentralapnea(k+1)==1
                Ecentralapnea(k)=1;
            end
            
        end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Clear Events which have been marked as central events
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    for k=1:length(BB_i_start)
        if (Ecentralhypop(k)||Ecentralapnea(k))&&removecentralhypops
            E1(k)=1;
%         elseif (Ecentralapnea(k))&&removecentralaps
%             E1(k)=1;
        end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

