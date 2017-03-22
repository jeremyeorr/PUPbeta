% This function implements the arousal shift and notation logic for the
% LGfromFlow_plot_V12_4 software. Depending on the option, different
% arousal markings are made. 


function [ AR1, numArousalBreaths, E_Terminate] = ArousalShift_V14_0( ArousalShift, AR1, E1, E_Terminate,BB_i_start ,numArousalBreaths)
       donotdeletearousalsinevents=1;
       % this stops the removal of arousals in events, implemented only for pi/4 and pi/16


       if (ArousalShift==0) 
            numArousalBreaths=1;
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                end
            end
        end

        %******************************************************************
%         This option sets all arousals as single breath arousals, and
%         shifts event termination arousals to the second breath following
%         the end of the events. Other arousals are shifted right by a
%         single breath
        %******************************************************************
        if (ArousalShift==pi/2) % 
            numArousalBreaths=2;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and set all arousals to a single breath,
            %starting after the first breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-1
                        bn=bn+1;
                    end
                    AR1(k:bn-1)=0;
                    AR1(k+1)=1;
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            %Check if there is an arousal in the period, and if so, define
            %an arousal on the second breath after the end of the event. 
            for k=1:(length(BB_i_start)-3)

                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                    if sum(AR1(k:k+3))>0
                        AR1(k:k+3)=0;
                        AR1(k+2)=1;
                    end

                end
            end
        end
        %******************************************************************
        %******************************************************************
        
        %******************************************************************
        % This option sets all arousals as single breath arousals, and
        % shifts event termination arousals to the second breath following
        % the end of the events. Other arousals are left at the original
        % location.
        %******************************************************************
        if (ArousalShift==pi/2.1) % 
            numArousalBreaths=2;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and set all arousals to a single breath,
            %starting after the first breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-1
                        bn=bn+1;
                    end
                    AR1(k:bn-1)=0;
                    AR1(k)=1;
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            %Check if there is an arousal in the period, and if so, define
            %an arousal on the second breath after the end of the event. 
            for k=1:(length(BB_i_start)-3)

                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                    if sum(AR1(k:k+3))>0
                        AR1(k:k+3)=0;
                        AR1(k+2)=1;
                    end

                end
            end
        end
        %******************************************************************
        %******************************************************************
        
        %******************************************************************
        % Leave arousals marked for their current duration, and apply a
        % VRA to each individual arousal marking.
        %******************************************************************
        if (ArousalShift==1) %was pi/4
            numArousalBreaths=1;
            
            if donotdeletearousalsinevents
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end
            end
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                end
            end
        end
        %******************************************************************
        %******************************************************************
        
        %******************************************************************
        % This option sets a 2 breath arousal parameter, where the second
        % breath effect is applied only to arousals which are marked for
        % more than 2 breaths. At event terminations, these are applied to
        % the second and third breaths following the end of an event...
        %******************************************************************
        if (ArousalShift==pi/8) % 
            numArousalBreaths=2.5;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and mark all arousals with the notation of 1 for the first breath of arousal, and 2 for all subsequent breaths. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-1
                        bn=bn+1;
                    end
                    AR1(k)=1; % Set first breath to arousal parameter 1.
                    AR1(k+1:bn-1)=2; % Set second, and any subsequent breaths to parameter 2. 
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            %Check if there is an arousal in the period, and if so, define
            %an arousal on the second breath after the end of the event. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                    if sum(AR1(k:k+3))>0
                        AR1(k:k+3)=0;
                        AR1(k+2)=1;
                        AR1(k+3)=2;
                    end
                end
            end
        end
        %******************************************************************
        %******************************************************************
        
        %******************************************************************
        % This option sets a 1 breath arousal parameter, where gamma1 is
        % applied to only first breath of any marked arousal.
        %******************************************************************
        if (ArousalShift==pi/14) % 
            numArousalBreaths=1;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and mark all arousals with the notation of 1 for any marked breath of the arousal, and 2 for a single "hangover" breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-2
                        bn=bn+1;
                    end
                    AR1(k)=1;
                    AR1(k+1:bn-1)=0; % Set any arousal marked breath to gamma 1
                    AR1(bn)=0; % Set breath following arousal to gamma2
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                end
            end
            
        end
        %******************************************************************
        %******************************************************************
        
        
        %******************************************************************
        % This option sets a 2 breath arousal parameter, where gamma1 is
        % applied to the first breath of any marked arousal, and gamma2 is
        % applied to any following marked breaths. No special adjustment at
        % event terminations.
        %******************************************************************
        if (ArousalShift==pi/15) % 
            numArousalBreaths=2.5;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and mark all arousals with the notation of 1 for any marked breath of the arousal, and 2 for a single "hangover" breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-2
                        bn=bn+1;
                    end
                    AR1(k)=1;
                    AR1(k+1:bn-1)=2; % Set any arousal marked breath to gamma 1
                    AR1(bn)=0; % Set breath following arousal to gamma2
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                end
            end
            
        end
        %******************************************************************
        %******************************************************************
        
        
        %******************************************************************
        % This option sets a 2 breath arousal parameter, where gamma1 is
        % applied to any marked arousal, and gamma2 is applied to the
        % breath following arousal (i.e. hangover breath). No special
        % adjustment at event terminations. 
        %******************************************************************
        if (ArousalShift==pi/16) % 
            numArousalBreaths=2.5;
            if donotdeletearousalsinevents
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            end
            %Scan through, and mark all arousals with the notation of 1 for any marked breath of the arousal, and 2 for a single "hangover" breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-2
                        bn=bn+1;
                    end
                    AR1(k:bn-1)=1; % Set any arousal marked breath to gamma 1
                    AR1(bn)=2; % Set breath following arousal to gamma2
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                end
            end
            
        end
        %******************************************************************
        %******************************************************************
        
        
        %******************************************************************
        % This option sets a 2 breath arousal parameter, where gamma1 is
        % applied to any marked arousal, and gamma2 is applied to the
        % breath following arousal (i.e. hangover breath). At event
        % terminations, these are applied to the second and third breaths
        % following the end of an event...
        %******************************************************************
        if (ArousalShift==pi/17) % 
            numArousalBreaths=2.5;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and mark all arousals with the notation of 1 for any marked breath of the arousal, and 2 for a single "hangover" breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-2
                        bn=bn+1;
                    end
                    AR1(k:bn-1)=1; % Set any arousal marked breath to gamma 1
                    AR1(bn)=2; % Set breath following arousal to gamma2
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            %Check if there is an arousal in the period, and if so, define
            %an arousal on the second breath after the end of the event. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                    if sum(AR1(k:k+3))>0
                        AR1(k:k+3)=0;
                        AR1(k+2)=1;
                        AR1(k+3)=2;
                    end
                end
            end
        end
        %******************************************************************
        %******************************************************************
        
        %******************************************************************
        % This option sets a 2 breath multiplicative arousal parameter,
        % where the second breath effect is applied only to arousals which
        % are marked for more than 2 breaths. At event terminations, these
        % are applied to the second and third breaths following the end of
        % an event...
        %******************************************************************
        if (ArousalShift==3) % 
            numArousalBreaths=3.5;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and mark all arousals with the notation of 1 for the first breath of arousal, and 2 for all subsequent breaths. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-1
                        bn=bn+1;
                    end
                    AR1(k)=1; % Set first breath to arousal parameter 1.
                    AR1(k+1:bn-1)=2; % Set second, and any subsequent breaths to parameter 2. 
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            %Check if there is an arousal in the period, and if so, define
            %an arousal on the second breath after the end of the event. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                    if sum(AR1(k:k+3))>0
                        AR1(k:k+3)=0;
                        AR1(k+2)=1;
                        AR1(k+3)=2;
                    end
                end
            end
        end
        
        %******************************************************************
        % This option sets a 2 breath multiplicative arousal parameter, where gamma1 is
        % applied to any marked arousal, and gamma2 is applied to the
        % breath following arousal (i.e. hangover breath). No special
        % adjustment at event terminations. 
        %******************************************************************
        if (ArousalShift==3.1) % 
            numArousalBreaths=3.5;
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and mark all arousals with the notation of 1 for any marked breath of the arousal, and 2 for a single "hangover" breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(AR1)-2
                        bn=bn+1;
                    end
                    AR1(k:bn-1)=1; % Set any arousal marked breath to gamma 1
                    AR1(bn)=2; % Set breath following arousal to gamma2
                    k=bn;
                end
                k=k+1;
            end
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)
                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
                    E_Terminate(k+1)=1;
                end
            end
            
        end
        %******************************************************************
        %******************************************************************
        
        if (ArousalShift==exp(pi)) % 
            AR1(1:end)=0; % Delete all arousals...
            
            %If current breath is an event, but the next 3 are not events,
            %define an end of breath termination. 
            for k=1:(length(BB_i_start)-3)

                if (E1(k)==0)&&E1(k+1)==1&&E1(k+2)==1&&E1(k+3)==1 
  
                    
                    if E_recover(k+1)==1
                        AR1(k+1)=1;
                    else
                        AR1(k+2)=1;
                    end
                end
            end
        end
        
        if (ArousalShift==pi)
            
            %If an arousal marked, and current and next breath are events, delete arousal..
            for k=1:(length(BB_i_start)-1)
                %if AR1(k)&&E_recover(k)==0&&E_recover(k+1)==0 % 
                if AR1(k)&&E1(k)==0&&E1(k+1)==0 %    
                    AR1(k)=0;
                end

            end
            %Scan through, and set all arousals to a single breath,
            %starting after the first breath. 
            k=1;
            while k<=(length(BB_i_start))
                if AR1(k)
                    bn=k;
                    while AR1(bn)&&bn<=length(BB_i_start)
                        bn=bn+1;
                    end
                    AR1(k:bn-1)=0;
                    AR1(k+1)=1;
                    k=bn;
                end
                k=k+1;
            end
        end
        
end

