function downloadRWS_URL(years,months,station,region,saveFolder)
% Save the RWS html files from UNIWYO for selected time and location



for Y=years
    for M=months
        while 1
            eom= num2str(eomday(Y, M));
            url=['http://weather.uwyo.edu/cgi-bin/sounding?region=' region '&TYPE=TEXT%3ALIST&YEAR=' num2str(Y) '&MONTH=' num2str(M) '&FROM=0100&TO=' eom '12&STNM=' station];
            filename=[saveFolder,'\RWS_' station '_' datestr(datenum(Y,M,1),'yyyy_mm') '.htm'];
            urlwrite(url,filename);
            % Control of file dimension!
            DR=dir(filename);
            if DR.bytes>500
                disp('Ok...')
                break%Continue only if dimension exceed 550 bytes. In case not, restart the same request
            else
                disp('Server too busy! Try again in 100 s')
                pause(100)
            end
        end
%     pause
pause(1)
    end
end