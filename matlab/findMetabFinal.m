function [outputMetab] = findMetabFinal(NMRDAT, NMRPAR, varargin)
%% AUTOMATED ANALYSIS OF ULTRA-HIGH RESOLUTION 2D-HSQC NMR SPECTRA
%
% The software automatically localizes a desired metabolite in your 2D-HSQC
% NMR spectra with an associated confidence level.
% It also estimates the contribution of each costituing part of a multiplet to the resulting pattern.
%
% FIRSTLY, launch the graphical user interface MetaboLab by the command sm
% and load your spectrum
%
%==========================================================================
% INPUT ARGUMENTS:
% NMRDAT                (struct)
% NMRPAR                (struct)
%==========================================================================
% OPTIONAL PARAMETERS:
% 'metabolites'         (Cell array) metabolites{:,1} the names of the metabolites (string)
%                                    metabolites{:,2} vector containing the spin numbers
% 'debug'               (string) Plot the output of FastICA algorithm and
%                       of the final regression (experimental data and
%                       estimated profile)
%                       Default is 'off'.
% 'report'              (string) Display a report (.html) containing
%                       figuers and the output of findMetab()
%                       Default is 'off'.
% 'plotFig'             (string) Plot a few figures during analysis
%                       Default is 'off'.
%==========================================================================
% OUTPUT ARGUMENTS:
% outputMetab           (struct)
%
%======================================================================
% EXAMPLES
% [output_find] = findMetab (dataHSQC, 'metabolite', yourCellArray ,'debug','on','report','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input requirement and optional parameters
if nargin <= 1
    error ('You must supply the HSQC data as input argument.');
end
sss = NMRPAR.CURSET(1);
eee = NMRPAR.CURSET(2);

% Set default values for optional parameters
metabolites = cell(0,0);
debug  = 'off';
report = 'off';
plotFig   = 'off';

maxWidth1H           = 0.1;
maxWidth13C          = 1.0;
minCorr              = 0.3;
c13threshold         = 5.0;
maxRange             = 1;
nReps                = 15;
scaleR2              = 1;
smooth               = 'off';
outputDir            = '.';
outputName           = 'autoHsqcMaReport';
xLibOffset           = 0.0;
yLibOffset           = 0.0;
rangeC               = 0.8;
nPlots               = 0;
r2value              = 1.5;
r2scaled             = 1.5;
index_good           = 0.0;

if (rem(length(varargin),2)==1)
    error('Provide the optional parameters by pair, name and option.');
else
    for i=1:2:(length(varargin)-1)
        % change the value of parameter
        switch lower (varargin{i})
            case 'metabolites'
                metabolites = varargin{i+1};
            case 'debug'
                debug = lower (varargin{i+1});
            case 'report'
                report = lower (varargin{i+1});
            case 'plotfig'
                plotFig = lower (varargin{i+1});
            case 'maxwidth1h'
                maxWidth1H = varargin{i+1};
            case 'maxwidth13c'
                maxWidth13C = varargin{i+1};
            case 'mincorr'
                minCorr = varargin{i+1};
            case 'c13threshold'
                c13threshold = varargin{i+1};
            case 'maxrange'
                maxRange = varargin{i+1};
            case 'nreps'
                nReps = varargin{i+1};
            case 'smooth'
                smooth = varargin{i+1};
            case 'outputdir'
                outputDir = varargin{i+1};
            case 'outputname'
                outputName = varargin{i+1};
            case 'xliboffset'
                xLibOffset = varargin{i+1};
            case 'yliboffset'
                yLibOffset = varargin{i+1};
            case 'rangec'
                rangeC = varargin{i+1};
            case 'nplots'
                nPlots = varargin{i+1};
            case 'r2value'
                r2value = varargin{i+1};
            case 'r2scaled'
                r2scaled = varargin{i+1};
            otherwise
                varargin{i}
                error ('Unknown optional parameter name.');
        end
    end
end

cnst        = find_var('$CNST=');
maxWidth13C = maxWidth13C*round(0.8*cnst(19)); %sqrt(cnst(19))); %floor(sqrt(cnst(19)));
maxRange    = maxRange + floor(cnst(19)/2);
dist1H      = 0.0;
dist13C     = 0.0;

switch debug
    case 'on'
        debugOn = 1;
    case 'off'
        debugOn = 0;
    otherwise
        error ('Unknown optional parameter value.');
end

switch report
    case 'on'
        reportOn = 1;
    case 'off'
        reportOn = 0;
    otherwise
        error ('Unknown optional parameter value.');
end

switch plotFig
    case 'on'
        plotFigOn = 1;
    case 'off'
        plotFigOn = 0;
    otherwise
        error ('Unknown optional parameter value.');
end

%% Set some variables
data = NMRDAT(sss,eee).MAT;
bool_date = 1;

%% Define the search area
size_list = size(metabolites,1);
for i = 1:size_list
    
    met  = str2func(metabolites{i,1});
    met  = met();
    
    if ismember(0,metabolites{i,2})
        for j =1:met.nspins
            spin(j,1) = struct('spin',struct('H1',zeros(1,2),'C13delta',[],'p_stim',[],'confidence',zeros(1,1)));
        end
    else
        for j =1:length(metabolites{i,2})
            spin(j,1) = struct('spin',struct('H1',zeros(1,2),'C13delta',[],'p_stim',[],'confidence',zeros(1,1)));
        end
    end
    outputMetab.(metabolites{i,1}) = spin;
    clear spin
end

%% Main loop
for iter_list = 1:size_list
    
    disp2(['Processing  Set: ' num2str(sss) ', Exp: ' num2str(eee) ', Metabolite: ' metabolites{iter_list,1}])
    if(NMRPAR.cmd==1)
        read_text_file
        drawnow()
    end
    m = str2func(metabolites{iter_list,1});
    m = m();
    for spinCounter = 1:m.nspins
        m.spin(spinCounter).c13 = m.spin(spinCounter).c13 + yLibOffset;
        m.spin(spinCounter).h1  = m.spin(spinCounter).h1  + xLibOffset;
    end
    
    if ismember(0,metabolites{iter_list,2})
        m_s = 1:m.nspins;
    else
        m_s = metabolites{iter_list,2};
    end
    
    isfinal = 0;
    
    for iter_spin = 1:length(m_s)
        name_sim = strcat(m.name2,m.spin(m_s(iter_spin)).name(end));
        %% FasICA algorithm
        maxWidth1H_range = [1];
        maxWidth13C_range = [1];
        index_good_old = 0;
        index_good     = 0;
        for range_area1 = 1:length(maxWidth1H_range)
            for range_area2 = 1:length(maxWidth1H_range)
                if debugOn == 1
                    [~, ~, m_shift_C_ppm, m_shift_H_ppm, range_H, range_C,~] = setSearchArea(NMRDAT,m,m_s(iter_spin), name_sim, debugOn, maxWidth1H, maxWidth13C);
                else
                    [~, ~, m_shift_C_ppm, m_shift_H_ppm, range_H, range_C] = setSearchArea(NMRDAT,m,m_s(iter_spin), name_sim, debugOn, maxWidth1H, maxWidth13C);
                end
                % Filters experimental signal only if smooth is on. This does not probably help with the localization accuracy and may
                % lead to underestimated multiplet percentages values
                x_filtered = zeros(size(data,1),length(range_H));
                if strcmp(smooth,'on')
                    for i = 1: length(range_H)
                        x_filtered(:,i) = sgolayfilt(data(:,range_H(i)),3,5);
                    end
                else
                    for i = 1: length(range_H)
                        x_filtered(:,i) = data(:,range_H(i));
                    end
                end 
                %% MAIN LOOP
                % the ICA algorithm is repeated nReps times setting the number of ICA components to be estimated to 3,5,7,10
                nIcaVector = [10]; % changed from 10 08/11/2021
                N_vett     = repmat(nIcaVector,nReps,1);
                N_vett     = N_vett(:)'; 
                mat_ica = cell(length(nIcaVector)*nReps,2);
                for n = 1:length(N_vett)
                    % Estimate ICA components
                    [ica_tmp] = fastica(x_filtered(range_C(1):range_C(end),:)','approach','defl','numOfIC',N_vett(n),'verbose','off','displayMode','off','g', 'skew','mu',1);
                    [sizeN1, ~] = size(ica_tmp);
                    for itr_ica = 1:sizeN1
                        if max(ica_tmp(itr_ica,:)) <= 0
                            ica_tmp(itr_ica,:) = ica_tmp(itr_ica,:) - min(ica_tmp(itr_ica,:));
                        end
                    end
                    mat_ica{n,1} = zeros(size(ica_tmp,1),size(x_filtered,1));
                    
                    try
                        mat_ica{n,1}(:,range_C(1):range_C(end)) = ica_tmp;
                    catch
                        keyboard;
                    end
                    % For each ICA component perform template matching to experimental 1D Carbon signal of the 2D spectrum by crosscorrelation
                    if debugOn
                        figure
                        count = 1;
                    end
                    index_all = zeros(length(mat_ica{n,1}(:,1)),1);
                    corr_all = zeros(length(mat_ica{n,1}(:,1)),1);
                    for j = 1:length(mat_ica{n,1}(:,1))
                        if debugOn == 1
                            count = mod(count-1,2*N_vett(n)) + 1;
                        end
                        template = (mat_ica{n,1}(j,:))';
                        template = template/norm(template);
                        % localize ICA component along proton dimension
                        [index_match, corr_match] = Find_match(template,data,range_C,range_H);
                        if length(corr_match) > 1
                            keyboard;
                        end
                        index_all(j) = index_match;
                        corr_all(j) = corr_match;
                        if debugOn
                            % Plots debug figures
                            template_x = 1:length(template);
                            template_x    = points2ppm(template_x,NMRDAT(sss,eee).PROC(2).REF);
                            index_all_ppm = points2ppm(index_all(j),NMRDAT(sss,eee).PROC(1).REF);
                            subplot(length(mat_ica{n,1}(:,1)),2,count)
                            plot(template_x,template)
                            title(sprintf('ICA n°%d | predicted ^1H shift = %2f | correlation to data = %.2f',j,index_all_ppm, corr_all(j)));
                            set(gca, 'XDir','reverse')
                            data_x = 1:length(data(:,1));
                            data_x  = points2ppm(data_x,NMRDAT(sss,eee).PROC(2).REF);
                            subplot(length(mat_ica{n,1}(:,1)),2,count+1)
                            plot(data_x,data(:,index_all(j)))
                            title(sprintf('match to ICA n°%d',j));
                            set(gca, 'XDir','reverse')
                            count = count+2;
                        end
                    end
                    mat_ica{n,2} = corr_all; % vector containing the maximum correlation scores of each ICA
                    mat_ica{n,3} = index_all; % The 1D carbon spectrum at index_all(i) is where the maximum correlation corr_all(i) is obtained
                end
                chooseNvett = zeros(length(mat_ica),1);
                for rep_N = 1:size(mat_ica,1)
                    if isempty(mat_ica{rep_N,2})
                        mat_ica{rep_N,2} = 0;
                    end
                    chooseNvett(rep_N) = max(mat_ica{rep_N,2});
                    if chooseNvett(rep_N) > 1
                        warning('In this repetition, more than 1 ICA signal has the same maximum correlation value!');
                        if debugOn == 1
                            keyboard;
                        end
                    end
                end
                index_N = find(chooseNvett == max(chooseNvett));
                if length(index_N) > 1 && debugOn == 1
                    warning('More than 1 ICA signal has the same maximum correlation value across repetitions!');
                    keyboard;
                end
                % Set of ICA signals of the repetition which contains at least one
                % ICA components with highest correlation 
                icasig = mat_ica{index_N(1),1}; 
                corr_all = mat_ica{index_N(1),2};
                index_all = mat_ica{index_N(1),3};

                mC13Shifts = zeros(1,length(icasig(:,1)));
                for ind_ic = 1:length(icasig(:,1))
                    [index_all(ind_ic), mC13Shifts(ind_ic)] = find_h_shift(icasig(ind_ic,:), index_all(ind_ic), data);
                end

                %% Multiplets simulation & fitting
                [s, c13shifts] = sim_multiplet(name_sim);

                gamma_adjust = 20.0*(42.577/10.7084);
                set_comb = 0;
                y_stim_NLS = zeros(length(icasig(1,:)),length(icasig(:,1)));
                R2_fitting = zeros(length(icasig(:,1)),1);
                peak_dec = zeros(length(icasig(:,1)),1);
                for i = 1:size(icasig,1)
                    shift_x = points2ppm(index_all(i),NMRDAT(sss,eee).PROC(1).REF);
                    y = icasig(i,:)';
                    y_temp = y/norm(y);
                    normS = norm(s(:,1));
                    s_1 = zeros(size(s));
                    for ii = 1:length(s(1,:))
                        s_1(:,ii) = (s(:,ii)/normS)*(max(y_temp)-min(y_temp));
                    end
                    [s_align] = alignment_multiplet((y_temp),(y),(s_1),range_C,set_comb,isfinal, minCorr, maxRange);
                    yy = data(:,index_all(i)); 
                    % Set to zero values outside the carbon and proton range
                    yy(1:range_C(1)-1) = zeros(length(yy(1:range_C(1)-1)),1);
                    yy(range_C(end)+1:end) = zeros(length(yy(range_C(end)+1:end)),1);
                    [~,y_NLS,R2] = estimate_LS(s_align,yy/norm(yy));
               
                    R2_fitting(i,1) = R2;
                    y_stim_NLS(:,i) = y_NLS;
                    peaks = peakdet(y_temp,0.5*max(abs(y_temp)));
                    peaks2 = ismember(peaks(:,1), range_C(1):range_C(end));
                    peak_inrange = ismember(round(mean(peaks(peaks2,1))),range_C(1) : range_C(end));
                    peaks_C_ppm(i,1) = points2ppm(round(mean(peaks(peaks2,1))),NMRDAT(sss,eee).PROC(2).REF);
                    dist2D(i,1)  = sqrt((shift_x - m_shift_H_ppm).^2 + ((peaks_C_ppm(i,1) - m_shift_C_ppm)/gamma_adjust).^2);
                    dist1H(i,1)  = abs(shift_x - m_shift_H_ppm);
                    dist13C(i,1) = abs(peaks_C_ppm(i,1) - m_shift_H_ppm);

                    if find(peak_inrange == 0) >= 1
                        peak_dec(i,1) = 0;
                    else
                        peak_dec(i,1) = length(peakdet(icasig(i,:),0.5*max(icasig(i,:)))); % how many peaks per ica exper
                    end

                end

              
                corr_high = find(corr_all(:) > minCorr);
                if isempty(corr_high)
                    for j = 1:length(corr_all)
                        peak_less(j,1) = length(peakdet(icasig(j,:),0.5*max(icasig(j,:))));
                    end
                    corr_high = find(min(peak_less) == peak_less);
                end

                if(length(corr_high)>length(peak_dec))
                    corr_high = corr_high(corr_high<=length(peak_dec));
                end
                for j = 1:length(corr_high)
                    try
                        peak_new(j,1) = peak_dec(corr_high(j),1); %num of peaks for high corr
                    catch
                        if(length(corr_high)>length(peak_dec))
                            corr_high = corr_high(1:length(peak_dec))
                        end
                        corr_high(corr_high>length(corr_high)) = length(corr_high);
                        peak_new(j,1) = peak_dec(corr_high(j),1);
                    end
                    index_new(j,1) = corr_high(j);
                end

                count_peak = 1;
                R2_fitting_high = [];
              
                for i = 1:length(corr_high)
                    R2_fitting_high(count_peak,:) = [R2_fitting(corr_high(i)) corr_high(i) R2_fitting(corr_high(i))*(1/dist2D(corr_high(i)).^2)]; %*max(abs(NMRDAT(cs,ce).MAT(range_C,index_all(corr_high(i)))))];
                    count_peak = count_peak +1;
                end
                if isempty(R2_fitting_high)
                    R2_fitting_high = [];
                end
                index_choose1 = R2_fitting_high(find(max(R2_fitting_high(:,3)) == R2_fitting_high(:,3)),2);
                %% Find the best 1H shift
                hill_climbing_on = 1;
                if hill_climbing_on 
                    [~,mIdx] = max(icasig(index_choose1,:));
                    h1sig     = data(mIdx,:);
                    diffSig   = -1;
                    sigIdx    = index_all(index_choose1);
                    idxRange = 1:15;
                    while(diffSig<0)
                        diffSig1 = h1sig(sigIdx) - h1sig(sigIdx-idxRange);
                        diffSig2 = h1sig(sigIdx) - h1sig(sigIdx+idxRange);
                        if(min(diffSig1)<0)
                            minIdx  = find(diffSig1==min(diffSig1)); %sigIdx - 1;
                            diffSig = diffSig1(minIdx); % = diffSig1;
                            sigIdx  = sigIdx - minIdx;
                        else
                            if(min(diffSig2)<0)
                                maxIdx  = find(diffSig2==min(diffSig2)); %sigIdx + 1;
                                diffSig = diffSig2(maxIdx); % = diffSig2;
                                sigIdx  = sigIdx + maxIdx;
                            else
                                diffSig = 1;
                            end
                        end
                    end
                    index_all(index_choose1) = sigIdx;
                end
                %% Final regression
                isfinal = 1;
                y_final = [];
                index_choose = index_choose1(1);
                y_final(:,1) = data(:, index_all(index_choose));
                y_final_filt = sgolayfilt(y_final,3,5);

                final_H = points2ppm(index_all(index_choose),NMRDAT(sss,eee).PROC(1).REF);
                y_final(1:range_C(1)-1) = zeros(length(y_final(1:range_C(1)-1)),1);
                y_final(range_C(end)+1:end) = zeros(length(y_final(range_C(end)+1:end)),1);

                normS = norm(s(:,1));
                for i = 1:length(s(1,:))
                    s(:,i) = (s(:,i)/normS)*(max(y_final)-min(y_final));
                end
                set_comb =1;

                y_final2 = abs(hft(y_final));
                peaks = peakdet(y_final2,0.5*max(y_final2));
                [Outliers] = find_outliers (peaks(:,1));
                if ~isempty(Outliers)

                    for j = 1:length(Outliers)
                        if ismember(Outliers(j),range_C)
                            y_final(Outliers(j)-1:Outliers(j)+1) = zeros(length(Outliers(j)-1:Outliers(j)+1),1);
                        end
                    end
                end

                isfinal = 1;

                index_choose1 = unique(index_choose1);
                if(~isvector((icasig(index_choose1,:)/norm(icasig(index_choose1,:)))'))
                    index_choose1 = index_choose1(1);
                end
                if(nPlots > 0)
                    
                    figHandle = figure;
                    n = index_choose1;
                    count = 1;
                    vett_plots =  1:nPlots;
                    if ~ismember(n,vett_plots)
                        vett_plots(1) = n;
                    end
                    for j = 1:length(vett_plots) 
                        % Plots debug figures
                        template_x = 1:length(template);
                        template_x    = points2ppm(template_x,NMRDAT(sss,eee).PROC(2).REF);
                        index_all_ppm = points2ppm(index_all(vett_plots(j)),NMRDAT(sss,eee).PROC(1).REF);
                      
                        subplot(nPlots,2,count)
                        plot(template_x,icasig(vett_plots(j),:)); %template)
                        set(gca, 'XLim', [template_x(range_C(end)), template_x(range_C(1))])
                        title(sprintf('ICA n°%d | predicted ^1H shift = %2f | correlation to data = %.2f',vett_plots(j),index_all_ppm, corr_all(vett_plots(j))));
                        set(gca, 'XDir','reverse')
                        data_x = 1:length(data(:,1));
                        data_x  = points2ppm(data_x,NMRDAT(sss,eee).PROC(2).REF);
                        subplot(nPlots,2,count+1)
                        plot(data_x,data(:,index_all(vett_plots(j))))
                        set(gca, 'XLim', [template_x(range_C(end)), template_x(range_C(1))])
                        title(sprintf('match to ICA n°%d',vett_plots(j)));
                        set(gca, 'XDir','reverse')
                        count = count+2;
                    end
                    set(figHandle,'PaperPositionMode','auto');         
                    set(figHandle,'PaperOrientation','landscape');
                    print([outputDir{1} filesep m.name2 num2str(iter_spin) '.pdf'], '-dpdf');
                    close(figHandle)
                end
                [s_align_final,~,~,D,D_best] = alignment_multiplet(y_final/norm(y_final),y_final,s,range_C,set_comb,isfinal,minCorr,maxRange);
               
                if(sum(s_align_final(:,end))==0)
                    if(sum(y_final/norm(y_final))<0)
                        if(max(y_final)==0)
                            y_final(round(length(y_final)/2)) = 1;
                            y_final(round(length(y_final)/2)+1) = 0.5;
                            y_final(round(length(y_final)/2)-1) = 0.5;
                        end
                        y_final(y_final<max(y_final)) = 0;
                    end
                    [s_align_final,~,~,D,D_best] = alignment_multiplet(y_final/norm(y_final),y_final,s,range_C,set_comb,isfinal,minCorr,maxRange);
                end
                for dd = 1:length(D_best)
                    D_best(dd) = D_best(dd) +D;
                end
             
                if(max(s_align_final(:,end))==0)
                    idx2 = find(sum(s_align_final)>0);
                    idx2 = idx2(length(idx2));
                    index_mult  = peakdet(s_align_final(:,idx2),0.5*max(s_align_final(:,idx2)));
                else
                    index_mult  = peakdet(s_align_final(:,end),0.5*max(s_align_final(:,end)));
                end
                ss          = NMRPAR.CURSET(1);
                ee          = NMRPAR.CURSET(2);
                rangePPM    = maxWidth13C; 
                rref        = NMRDAT(ss,ee).PROC(2).REF;
                rangePoints = ppm2points([0, rangePPM],rref);
                rangePoints = round(abs(diff(rangePoints)));

                index_start = index_mult(1,1) - rangePoints;
                index_end = index_mult(end,1) + rangePoints;
                y_final_spc = y_final/norm(y_final);
            
                if(~isnan(y_final_spc))
                    if(0)
                        disp('////////////////////////////////////////////////////////////////////////////////////')
                        [s_align_final2_0,~,~,D2_0,D_best2_0] = alignment_multiplet(y_final_spc,y_final,s_align_final,range_C,set_comb,isfinal,minCorr,0.0);
                        for dd = 1:length(D_best2_0)
                            D_best2_0(dd) = D_best(dd) + D_best2_0(dd) +D2_0;
                        end
                        [p_LS_0, y_stim_final_0,index_good_0,index_good_component_0] = estimate_LS(s_align_final2_0,y_final);
                        [s_align_final2_mr,~,~,D2_mr,D_best2_mr] = alignment_multiplet(y_final_spc,y_final,s_align_final,range_C,set_comb,isfinal,minCorr,maxRange);
                        for dd = 1:length(D_best2_mr)
                            D_best2_mr(dd) = D_best(dd) + D_best2_mr(dd) +D2_mr;
                        end
                        [p_LS_mr, y_stim_final_mr,index_good_mr,index_good_component_mr] = estimate_LS(s_align_final2_mr,y_final);
                        if(sum(index_good_component_0) > sum(index_good_component_mr))
                            s_align_final2       = s_align_final2_0;
                            D2                   = D2_0;
                            D_best2              = D_best2_0;
                            p_LS                 = p_LS_0;
                            y_stim_final         = y_stim_final_0;
                            index_good           = index_good_0;
                            index_good_component = index_good_component_0;
                        else
                            s_align_final2       = s_align_final2_mr;
                            D2                   = D2_mr;
                            D_best2              = D_best2_mr;
                            p_LS                 = p_LS_mr;
                            y_stim_final         = y_stim_final_mr;
                            index_good           = index_good_mr;
                            index_good_component = index_good_component_mr;
                        end
                    else
                       
                        [s_align_final2,~,~,D2,D_best2] = alignment_multiplet(y_final_spc,y_final,s_align_final,range_C,set_comb,isfinal,minCorr,maxRange);
                        for dd = 1:length(D_best2)
                            D_best2(dd) = D_best(dd) + D_best2(dd) +D2;
                        end
                        [p_LS, y_stim_final,index_good,index_good_component] = estimate_LS(s_align_final2,y_final);
                  
                    end

                   
                    % Penalise singlet peaks
                    maxDist = sqrt(maxWidth1H.^2 + (maxWidth13C/gamma_adjust).^2);
                    weightedDistance = 1;
                    if dist2D(index_choose1) > maxDist
                        weightedDistance = weightedDistance * (1 - (dist2D(index_choose1)/maxDist - 1));
                    else
                        weightedDistance = weightedDistance * 1;
                    end
                    if(p_LS(1) == 1)
                        index_good = index_good - 0.5;
                    end
                    if(0.5*max(y_final(range_C))<=0)
                        y_final(round(mean(range_C))) = 1;
                        y_final(round(mean(range_C))+1) = 0.5;
                        y_final(round(mean(range_C))-1) = 0.5;
                    end
                    peaks_R2 = peakdet(y_final(range_C),0.5*max(y_final(range_C)));

                    if size(peaks_R2,1)>5
                        index_good = 0;
                    end
                    p_LS = p_LS/norm(p_LS);
                    peaks_final = [];
                    peaks_final_stim = [];
                    peaks_final_ppm = [];
                    peaks_final_stim_ppm = [];
                    peaks_final = peakdet(y_final,0.5*max(y_final));
                    if(max(y_stim_final) <= 0)
                        peaks_final_stim = zeros(size(peaks_final));
                    else
                        peaks_final_stim = peakdet(y_stim_final,0.5*max(y_stim_final));
                    end

                    for k= 1:length(peaks_final(:,1))
                        peaks_final_ppm(k,1) = points2ppm(peaks_final(k,1),NMRDAT(sss,eee).PROC(2).REF);
                    end
                    for k= 1:length(peaks_final_stim(:,1))
                        peaks_final_stim_ppm(k,1) = points2ppm(peaks_final_stim(k,1),NMRDAT(sss,eee).PROC(2).REF);
                    end

                    y_stim_final_x = range_C(1):range_C(end);
                    for k= 1:length(y_stim_final_x)
                        y_stim_final_x(k) = points2ppm(y_stim_final_x(k),NMRDAT(sss,eee).PROC(2).REF);
                    end

                    y_final_x = (range_C(1):range_C(end));
                    for k= 1:length(y_final_x)
                        y_final_x(k) = points2ppm(y_final_x(k),NMRDAT(sss,eee).PROC(2).REF);
                    end

                    if(strcmp(plotFig,'on'))
                        h = figure;
                        hold all
                        plot(y_final_x,y_final(range_C(1):range_C(end))/norm(y_final(range_C(1):range_C(end))))
                        plot(y_stim_final_x,y_stim_final(range_C(1):range_C(end))/norm(y_stim_final(range_C(1):range_C(end))))

                        legend('Data','estimation')
                        xlabel('^{13}C [ppm]');
                        set(gca, 'XDir','reverse')

                        title(sprintf('%s D = [1H = %d, 13C = %d] | P = [1H = %d] | R2 = %d',name_sim, m_shift_H_ppm,m_shift_C_ppm,final_H, index_good));
                    end

                    %
                    origin_dir = pwd;
                    if bool_date == 1
                        currDate = strrep(datestr(datetime), ':', '_');
                        currDate(isspace(currDate)) = [];
                        bool_date = 0;
                    end
                    if ~exist('Results','dir')
                        mkdir('Results',currDate)
                        my_result_fold = currDate;
                    end
                    if(0)
                        if exist('Results','dir') && bool_date == 0 && bool_folder == 0
                            cd Results
                            mkdir(currDate)

                            bool_folder = 1;
                            cd(currDate)
                            if(strcmp(plotFig,'on'))
                                saveas(h,sprintf('fig_%s.fig',name_sim))
                            end
                        else
                            cd Results
                            cd(currDate)
                            if(strcmp(plotFig,'on'))
                                saveas(h,sprintf('fig_%s.fig',name_sim))
                            end
                        end
                    end
                    cd(origin_dir);

                    if index_good >= index_good_old

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.H1          = [final_H index_all(index_choose)];

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.C13delta    = D_best2;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.p_stim      = p_LS/sum(p_LS);

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.confidence  = weightedDistance*index_good; % [15/11/2021 - adjust coefficient of determination according to distance from lilbrary peak]

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.confidence2 = index_good_component;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.C13         = c13shifts;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.simSpc      = y_stim_final;
                       
                        index_good_old = index_good;
                    else
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.H1          = [final_H index_all(index_choose)];
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.C13delta    = D_best2;
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.p_stim      = zeros(size(p_LS/sum(p_LS)));
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.confidence  = 0;
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.confidence2 = 0; 
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.C13         = c13shifts;
                    end
                else
                    if index_good >= index_good_old
                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.H1         = [0];

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.C13delta   = 0;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.p_stim     = 0;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.confidence = 0;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.confidence2 = 1;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.C13        = 80;

                        outputMetab.(metabolites{iter_list,1})(m_s(iter_spin)).spin.simSpc     = 0;
                        
                        index_good_old = index_good;
                    end
                end
            end
        end
    end
end
if(0)
    cd Results
    cd(currDate)
    save('outputMetab.mat','outputMetab')
    cd(origin_dir)
end
if reportOn == 1
    idxFind = strfind(outputDir,filesep);
    strLen  = length(outputDir);
    if(idxFind(length(idxFind))==strLen)
        outputDir = outputDir(1:strLen-1);
    end
    %% Generate report!
    hsqc_generate_report
end
return


