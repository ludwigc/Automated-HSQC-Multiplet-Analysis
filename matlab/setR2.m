function setR2(r2,fName,pName)

global NMRPAR
global NMRDAT

s = NMRPAR.CURSET(1);
e = NMRPAR.CURSET(2);

ref = NMRDAT(s,e).PROC(2).REF;


mlsys = read_mlsys(fName);

for k = 1:length(mlsys)
    mlsys(k).deltaR2 = r2;
end

make_mlsys2_hsqcma_ica(mlsys, ref, r2, [pName fName '.mlsys']);