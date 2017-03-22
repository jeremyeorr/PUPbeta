% This function is designed to fix clipping in a respiratroy flow signal-
% with a particular focus on the ACZ database. It uses a percentile base to
% identify fistly whther there is clipping, and then identifies the regions
% where clipping may occur. Interpolation is used to fill in the regions
% where clipping is present. To assist in correct interpolation, data is
% down sampled, and then interp1 function to avoid crazy high freqeuncy
% interpolation. Finally, the original signal is used outside clipping
% areas, and the upsampled interpolation used in the clipping resion. 
% 

function [x3,clipped] = FixClipping_V14_0(t1,x1,dt2)


    dt=(t1(end)-t1(1))/(length(t1)-1);
    %dt2/dt;
    tds=downsample(t1,round(dt2/dt));
    xds=downsample(x1,round(dt2/dt));

    delta=max(x1)-min(x1);
    maxthres=max(x1)-0.01*delta;
    minthres=min(x1)+0.01*delta;

    x2=xds;
    t2=tds;
    clippedds=(xds>=maxthres)|(xds<=minthres);
    x2(clippedds)=[];
    t2(clippedds)=[];
    
    if (length(t2)<=2)||(length(x2)<=2)
        clipped(1:length(x1))=1;
        x3=NaN;
    else
        x3=interp1(t2,x2,t1,'spline');
        clipped=interp1(tds,clippedds,t1,'linear');
        clipped(clipped>0)=1;
        x1(clipped>0)=x3(clipped>0);
    end
    
%     figure(1)
%     plot(t2,x2,'g',t1,[x3;clipped;x1]);    


end

