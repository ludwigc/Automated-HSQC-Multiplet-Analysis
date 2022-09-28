function [multiplets, c13shifts] = sim_multiplet(metName) 

if(nargin<1)
    metName = 'lactate3';
end
%if(nargin<2)
%    lb = 4;
%end

global NMRDAT
global NMRPAR

if(~isfield(NMRPAR,'closeTerminal'))
    NMRPAR.closeTerminal = 1;
end


c13shifts  = [];

s          = NMRPAR.CURSET(1);
e          = NMRPAR.CURSET(2);

origDir    = pwd;
mlsysDir = NMRPAR.mlsysDir;

cd(mlsysDir)
mlsys      = read_mlsys_ica(metName);
nSubSys    = length(mlsys);
[nPts,~]   = size(NMRDAT(s,e).MAT);
ref        = NMRDAT(s,e).PROC(2).REF;
td         = NMRDAT(s,e).ACQUS(2).TD/2;
angle      = NMRDAT(s,e).PROC(2).SSB;
int        = 1;
phase      = 0;

[val, valIdx] = find_var('$CNST=');
if(isempty(val))
    jresConst  = 1;
    jresNConst = 1;
else
    jresNConst = val(18);
    jresConst  = val(19);
    if(jresConst<0)
        jresConst = 1;
    end
    if(jresNConst<0)
        jresNConst = 1;
    end
end



if(isempty(angle))
    angle = 90;
    NMRDAT(s,e).PROC(2).SSB = 90;
end
if(ischar(NMRPAR.echoTime))
    echoTime = str2num(NMRPAR.echoTime);
else
    echoTime = NMRPAR.echoTime;
end
if(ischar(NMRPAR.echoTime2))
    echoTime2 = str2num(NMRPAR.echoTime2);
else
    echoTime2 = NMRPAR.echoTime2;
end

multiplets = zeros(nPts,nSubSys);

for k = 1:nSubSys
    for l = 1:nSubSys
        if(k==l)
            mlsys(l).mult = 1;
        else
            mlsys(l).mult = 0;
        end
    end
    make_mlsys_hsqc(mlsys, ref); %, lb);
    writebin_hsqcma('temp',ref,multiplets(:,k),td);
    [a,b]     = system([NMRPAR.pythonExec ' -V']);
    b         = strtrim(b);
    if(strfind(b,'2.'))
        simspcStr = replace(which('simspc.py'),' ','\ ');
    else
        simspcStr = replace(which('simspc_p3.py'),' ','\ ');
    end
    if(ismac)
        systemString = ['cd ' pwd() ';' NMRPAR.pythonExec ' ' simspcStr ' ' num2str(int) ' ' num2str(phase) ' ' num2str(angle) ' ' num2str(echoTime/jresConst) ' ' num2str(echoTime2/jresNConst) ' ' num2str(jresConst) ' ' num2str(jresNConst) '; exit'];
        [a,b] = system(systemString);
        if(0)
            terminalNotOpened = 1;
            [a,kk]  = system(['ps -U$USER|grep "Terminal"|grep -v grep']);
            if(~isempty(kk))
                terminalNotOpened = 0;
            end
            sysCmd2 = ['osascript -e ''tell Application "Terminal" to do script "' systemString '"''']; %'"; set newWindow to id of front window; tell window id newWindow, set index to 1, set visible to false, end tell'''];
            sysCmd2 = replace(sysCmd2,'\ ','\\ ');
            [a,b]   = system(sysCmd2);
            kk      = 'abcde';
            while(~isempty(kk))
                pause(1)
                [a,kk] = system(['ps |grep "python"|grep -v grep']);
                pause(1)
            end
            if(NMRPAR.closeTerminal==1)
                sysCmd2 = ['osascript -e ''tell Application "Terminal" to close'''];
                system(sysCmd2);
                if(terminalNotOpened==1)
                    sysCmd2 = ['osascript -e ''tell Application "Terminal" to quit'''];
                    system(sysCmd2);
                end
            end
        end
    else
        systemString = ['system(''' NMRPAR.pythonExec ' ' simspcStr ' ' num2str(int) ' ' num2str(phase) ' ' num2str(angle) ' ' num2str(echoTime/jresConst) ' ' num2str(echoTime2/jresNConst) ' ' num2str(jresConst) ' ' num2str(jresNConst) ''');'];
        eval(systemString);
    end
    make_mlsys_hsqc(mlsys, ref, [metName '.mlsys']);
    multiplets(:,k) = readbin_temp;
end


c13shifts = zeros(1,length(mlsys));
for k = 1:length(mlsys)
    c13shifts(k) = mlsys(k).cs(1);
end

cd(origDir)


