%readxml from NSRR

filenamexml='S:\RESEARCH\STUDIES\3. DATA\Jeremy Orr\HeartBeat\heartbeat-baseline-700001-nsrr.xml'

S=xml2struct(filenamexml)

try
Nepochs = length(S.PSGAnnotation.SleepStages.SleepStage);
clear Epochs
Epochs=NaN*zeros(Nepochs,1);
for i=1:Nepochs
Epochs(i) = str2num(S.PSGAnnotation.SleepStages.SleepStage{1,i}.Text);
end %W=0,N1=1,N2=2,N3=3,R=5
figure(1); plot(Epochs);
end

Nevents = length(S.PSGAnnotation.ScoredEvents.ScoredEvent);
clear Events
for i=1:Nevents
Events{i}.name = S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.EventConcept.Text;
Events{i}.start = str2num(S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.Start.Text);
Events{i}.duration = str2num(S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.Duration.Text);
end %W=0,N1=1,N2=2,N3=3,R=5


Eventnames = {'Arousal (ARO RES)','Obstructive Apnea','Mixed Apnea','Hypopnea','Central Apnea'}

%find unique names
if 0
UniqueNames=[];
for i=1:Nevents
    unique=1;
    for n=1:length(UniqueNames)
        if strcmp(Events{i}.name.Text,UniqueNames{n})
            unique=0;
            break
        end
    end
    if unique
        UniqueNames{length(UniqueNames)+1}=Events{i}.name.Text
    end
end
end


Eventnames = {'Arousal (ARO RES)','Obstructive Apnea','Mixed Apnea','Hypopnea','Central Apnea'}

