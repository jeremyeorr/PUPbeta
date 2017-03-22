function getEDFchannelnames()
global ExportDataSpreadsheet handletext
[~,txt] = xlsread(ExportDataSpreadsheet,1,'E4:K10000');

if 0
    for n=1:size(txt,1)
        directory = char(txt(n,4));
        fname=char(txt(n,1));
        Label = EDFChannelLabels([directory '\' fname]);
        xlswrite(ExportDataSpreadsheet,Label,1,['U' num2str(3+n)]);
    end
else
    N=size(txt,1);
    Labels = cell(N,50); % Assuming max 50 channels
    for n=1:N%size(txt,1)
        directory = char(txt(n,4));
        fname=char(txt(n,1));
        try
            Label = EDFChannelLabels([directory '\' fname]);
            Labels(n,1:length(Label)) = Label;
            displaytext=['Reading n=' num2str(n) '/' num2str(N) ', ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
        catch me
            displaytext=['Reading n=' num2str(n) ', Error: ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
        end
    end
    displaytext=['Writing data to Xls'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    xlswrite(ExportDataSpreadsheet,Labels,1,['U' num2str(3+1)]);
    displaytext=['Complete'];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

function Label = EDFChannelLabels(fname)
fid = fopen(fname,'r');
fseek(fid,252,-1);
M  = sscanf(char(fread(fid,4,'char')'),'%d');
for m=1:M+12
    fseek(fid,(m-1)*16+256,-1);
    Label{m} = char(fread(fid,16,'char')'); %%%%%%%%%
end
fclose(fid); % Close file
