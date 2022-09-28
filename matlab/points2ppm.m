function cs=points2ppm(point,ref_ppm,ref_point,TD,SWH,SFO1)
%POINTS2PPM - Converts points to ppm if REF is available
%
%             call:  cs=points2ppm(point,ref_ppm,ref_point,TD,SWH,SFO1)
%
% See also: PPM2POINTS
  
% ULG - 30-05-01
% Copyright (c) U.L. G�ther  & C. Ludwig 2001

if(length(ref_ppm)>1)
    ref_point = ref_ppm(2);
    TD        = ref_ppm(3);
    SWH       = ref_ppm(4);
    SFO1      = ref_ppm(5);
    ref_ppm   = ref_ppm(1);
end


SW = SWH / SFO1;
D_ppm   = SW /(TD-1);
ppm_offs = ref_ppm - (TD-ref_point)*D_ppm;

% temp = linspace(SW,0,TD);	% from SW to 0 in TD steps
% axis = temp + ppm_offs;

% chemical shifts in ppm:
%ref_ppm
%ref_point
cs = ppm_offs + (TD-point)*D_ppm;

%if nargout~=0
%  ppm_axis = axis;
%end

return
