function [m_shift_C,m_shift_H,m_shift_C_ppm,m_shift_H_ppm,range_H,range_C,fig_c] = setSearchArea(NMRDAT,m,m_s,name_sim,debugOn,maxWidth1H,maxWidth13C)
%% Constrain the search space of the 2D NMR spectrum to an area around the theoretical position of the metabolite (1H,13C)
% INPUT: NMRDAT (NMR data); m (metabolite struct); m_s (metabolite's spin
% number); name_sim (metabolite multiplet name); debugOn (if == 1 plots a
% debug figure); maxWidth1H (max proton width in ppm); maxWidth13C (max carbon width in ppm);
% OUTPUT: m_shift_C (multiplet's central peak carbon shift in points); m_shift_H (multiplet's central peak proton shift in points);
% m_shift_C_ppm (multiplet's central peak carbon shift in ppm);
% m_shift_H_ppm (multiplet's central peak proton shift in ppm); range_H
% (proton range in terms of columns rage of the 2D NMR spectrum); range_C
% (carbon range in terms of raws rage of the 2D NMR spectrum); fig_c figure
% handle
%%
% get the shifts in the C and H dimension (points)
global NMRPAR
s = NMRPAR.CURSET(1);
e = NMRPAR.CURSET(2);
ref2 = NMRDAT(s,e).PROC(2).REF;
ref1 = NMRDAT(s,e).PROC(1).REF;
m_shift_C = ppm2points(m.spin(m_s).c13,ref2);
m_shift_H = ppm2points(m.spin(m_s).h1,ref1);
% get the shifts in the C and H dimension (ppm)
m_shift_C_ppm = m.spin(m_s).c13;
m_shift_H_ppm = m.spin(m_s).h1;
%find how many elements in the data set are in x ppm (depends on the resolution)
[elementsPerppm_H, elementsPerppm_C ] = resolution_ppm(maxWidth1H, maxWidth13C);
% Define search area
range_H = m_shift_H - elementsPerppm_H :1: m_shift_H + elementsPerppm_H;
range_C = m_shift_C - elementsPerppm_C :1: m_shift_C + elementsPerppm_C;
%% restrict data points within given maxWidth1H/13C
idxH = abs(points2ppm(range_H,ref1) - m.spin(m_s).h1)<maxWidth1H;
idxC = abs(points2ppm(range_C,ref2) - m.spin(m_s).c13)<maxWidth13C;
range_H = range_H(idxH);
range_C = range_C(idxC);
% Plot 3D
if debugOn == 1
    data = NMRDAT.MAT;
    range_C_ppm = zeros(length(range_C),1);
    for k= 1:length(range_C)
        range_C_ppm(k) = points2ppm(range_C(k),NMRDAT(s,e).PROC(2).REF);
    end
    range_H_ppm = zeros(length(range_H),1);
    for k= 1:length(range_H)
        range_H_ppm(k) = points2ppm(range_H(k),NMRDAT(s,e).PROC(1).REF);
    end
    [xx,yy] = meshgrid(range_H_ppm,range_C_ppm);
    i = im2double(data(range_C,range_H));
    fig_c = figure;
    contour(xx,yy,i)
    xlabel('^1H [ppm]');
    ylabel('^{13}C [ppm]');
    set(gca, 'XDir','reverse')
    set(gca, 'YDir','reverse')
    title(sprintf('%s, S/N = %4.2f',name_sim,max(max(abs(i)))/std(i(:,1))));
end
end