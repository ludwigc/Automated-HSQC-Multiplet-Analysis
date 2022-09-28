function [sigIdx, mC13Shift] = find_h_shift(ic, sigIdx, data)
    %% Find the best 1H shift
    [~,mIdx] = max(ic);
    global NMRDAT; 
    ref2=NMRDAT(1,1).PROC(2).REF; 
    mC13Shift = points2ppm(mIdx,ref2);
    h1sig     = data(mIdx,:);
    diffSig   = -1;
    while(diffSig<0)
        diffSig1 = h1sig(sigIdx) - h1sig(sigIdx-1);
        diffSig2 = h1sig(sigIdx) - h1sig(sigIdx+1);
        if(diffSig1<0)
            sigIdx  = sigIdx - 1;
            diffSig = diffSig1;
        else
            if(diffSig2<0)
                sigIdx = sigIdx + 1;
                diffSig = diffSig2;
            else
                diffSig = 1;
            end
        end
    end
end
       
        
        