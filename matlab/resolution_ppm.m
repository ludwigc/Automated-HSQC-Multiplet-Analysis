function [elementsPerHppm, elementsPerCppm, ppmPerPointH, ppmPerPointC] = resolution_ppm(maxWidth1H, maxWidth13C)

if(nargin<1)
    maxWidth1H = 0.15;
end
if(nargin<2)
    maxWidth13C = 0.7;
end
% 
global NMRDAT
global NMRPAR
cs = NMRPAR.CURSET;

max_width_H = maxWidth1H;
max_width_C = maxWidth13C;

cnst = find_var('$CNST='); %Bruker variable
cnst18 = cnst(19);
%
refH = NMRDAT(cs(1),cs(2)).PROC(1).REF;
refC = NMRDAT(cs(1),cs(2)).PROC(2).REF;

ppmPerPointH = abs(diff(points2ppm([1,2],refH)));
ppmPerPointC = abs(diff(points2ppm([1,2],refC)));

elementsPerHppm = round(max_width_H/ppmPerPointH);
elementsPerCppm = round((max_width_C*cnst18)/ppmPerPointC);  
end

