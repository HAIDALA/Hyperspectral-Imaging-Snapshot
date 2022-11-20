function res= load_data_from_spectrometer_file(path)

fnames=dir(path);

X = [];
for zz=1:numel(fnames)-2
    fname=fnames(zz+2).name;
    fprintf('%s, %d out of %d\n',fname,zz,numel(fnames)-2);
    %     fid =fopen(['C:\Users\Kinan\Desktop\LISIC\',fname]);
    fid =fopen([fname]);
    C=textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    x=[];
    for k=numel(C{1,1}):-1:1
        tmp=C{1,1}{k,1};
        
        if(contains(tmp,'Wavelength'))
            break;
        end
        tmp=tmp(5:end);
        tmp=str2num(tmp);
        x=[tmp;x];
    end
    X = [X,x];
    
end
fnames_new={};
for i=3:size(fnames)
    temp=strsplit(fnames(i).name,'.');
    fnames_new{i}=[char(temp{1,1}),char(temp{1,2}),char(temp{1,3})];
end

res = array2table(X,'VariableNames',fnames_new(3:end));
end

