function [dates,dat,headers] = gFileOpen(fileName,treatAsWtr)

if eq(nargin,1)
    treatAsWtr = false;
end
% author: Luke A Winslow March 2011
    fid = fopen(fileName);
    

    headers = fgetl(fid);
    format = '%s';
    for i=2:length(regexp(headers,'\t','split'));
        format = strcat(format,'%f');
    end
    d = textscan(fid,format,'delimiter','\t','treatAsEmpty',{'na','NA',...
        '#VALUE!','#NAME?'});
    fclose(fid);
    
    
    dat = horzcat(d{2:end});
    try dates = datenum(d{1},'yyyy-mm-dd HH:MM');
        
    catch mssg
        
        dates = NaN(length(d{1}),1);
        for i = 1:length(d{1})
            dates(i) = datenum(char(d{1}(i)),'yyyy-mm-dd HH:MM');
        end
    end
    
    if treatAsWtr
        heads = textscan(headers,'%s','Delimiter','\t');
        depth = NaN(1,length(heads{1})-1);
        for i = 1:length(heads{1})-1
            txtT = char(heads{1}(i+1));
            depth(i) = str2double(txtT(5:end));
        end
        headers = depth;  % NOW are numbers...
    end
    

end