% Running the full dissipation on the arc length resolved data.

 %%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 Drive_Letter = 'H';
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Constants

 Dj = 9.525/1000; % jet diameter (m)
 c = 2.12; %Compensated structure function constant
 
 
 % File stuff
 Lstr.rff_inst_pth = [Drive_Letter ':\JICF_Port\Drive_Data\Dissipation_RFF_ArcLength'];
 Lstr.diss_pth = [Drive_Letter ':\JICF_Port\Drive_Data\Dissipation_Multiple_Ways'];
 Lstr.nff_stat_pth = [Drive_Letter ':\JICF_Port\Drive_Data\MFStats_Update'];
 
 Lstr.rff_inst_nams = list_mats(Lstr.rff_inst_pth);
 Lstr.diss_nams = list_mats(Lstr.diss_pth);
 Lstr.nff_stat_nams = list_mats(Lstr.nff_stat_pth);
 
 Sstr.strPth = [Drive_Letter ':\JICF_Port\Drive_Data\Dissipation_Arc_Length_Data_Out\stFun'];
 Sstr.dissPth = [Drive_Letter ':\JICF_Port\Drive_Data\Dissipation_Arc_Length_Data_Out\diss_corr'];
 Sstr.probPth = [Drive_Letter ':\JICF_Port\Drive_Data\Dissipation_Arc_Length_Data_Out'];
 % find only diss_nams in SLRel
 Lstr.diss_nams = Lstr.diss_nams(contains(Lstr.diss_nams, 'SL'));
 
 Lstr.unders = cellfun(@(x) strfind(x, '_'), Lstr.diss_nams, 'UniformOutput', false);
 Lstr.names = cellfun(@(x,y) x(1:y(2)-1), Lstr.diss_nams, Lstr.unders,...
    'UniformOutput', false);

%%
fullTime = zeros(1,length(Lstr.names));
sumTime = 0;
for ii = 2:length(Lstr.names)  
    
    fullT = tic;
    
    nff_dat = load(char(strcat(Lstr.nff_stat_pth, filesep, Lstr.nff_stat_nams(...
        contains(Lstr.nff_stat_nams, Lstr.names{ii})))));
    inst_dat = load(char(strcat(Lstr.rff_inst_pth, filesep, Lstr.rff_inst_nams(...
        contains(Lstr.rff_inst_nams, Lstr.names{ii})))));
    diss_dat = load(char(strcat(Lstr.diss_pth, filesep, Lstr.diss_nams(...
        contains(Lstr.diss_nams, Lstr.names{ii})))));
    
    Vj = nff_dat.Uj;
    normFac = Vj^3./Dj;
    
    SCreal = inst_dat.SCmat .* Dj;
    %Pre-allocating the full structure function things
    stDx = zeros(200, size(inst_dat.XCmat, 2), size(inst_dat.Vr, 3));
    stDy = stDx;
    
    %Variables to pass the parallel loop 
    parNam = Lstr.names{ii};
    parTot = length(Lstr.names);
    parSCLen = size(inst_dat.Vr, 3);
    parVel = inst_dat.Vr;
    
    % Parallel loop runs the structure functions
    parfor jj = 1:size(inst_dat.Vr, 3)
       
        sfunT = tic;
        
        [stDx(:,:,jj), stDv(:,:,jj)] = PIV_Struc_FunY(SCreal, parVel(:,:,jj));
        time = toc(sfunT);
        
        IterTime = datestr(datenum(0,0,0,0,0,fullTime(ii)), 'HH:MM:SS');
        SumTime = datestr(datenum(0,0,0,0,0,sumTime), 'HH:MM:SS');
        fprintf('\n%s structure fun calc: \nIteration %d of %d took %3.2f seconds \n Run %d of %d, last run took %s\n Run time so far %s\n'...
            ,parNam, jj, parSCLen, time, ii, parTot, IterTime, SumTime);
        
    end
    
    % Takes the means of the structure fuctions 
    stFun.Mdx = nanmean(stDx, 3);
    stFun.Mdv = nanmean(stDv, 3);
    
    % Checks if there is anything in the structure functions, continues to
    % next itteration if no-data. Populates a 1 in the 'problem' structure
    % to indicate shits fucked up.
    if nnz(isnan(stFun.Mdv)) == numel(stFun.Mdv)
        problem.(Lstr.names{ii}) = 1;
        continue
    else
        problem.(Lstr.names{ii}) = 0;
    end
    
    %%%%%%%%%%% This is some weird normalization stuff %%%%%%%%%%%%%%%%%
    % If the normalization on the instantanious is fixed this should change
    % to normalizing the Mdv./Vj^2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stFun.Ndv = stFun.Mdv./Vj^4;
    stFun.Ndx = stFun.Mdx./Dj;
    
    % Compensated structure function
    stFun.Cdv = (stFun.Ndv./c).^(3/2).*(1./stFun.Ndx);
    
    % Find dissipation via structure functions with the Struct_Inflects function 
    [stFun.Sdiss, stFun.inStr] = Struc_Inflects(stFun.Cdv, stFun.Ndx);
    
    % Take the mean of the gradient methods at each column 
    diss_fields = {'D_Graham', 'D_Greg', 'D_Xu'};
    
    %Caclulate the modified structure function dissipation, this step is
    %done with the dimensional dissipation from the structure functions
    
    %%%%%%%%%%%% Note 3/20/20: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perhaps del value inside ModifiedStructureFun should be adjusted 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kk = 1:length(stFun.Sdiss)

        [diss_corr.Rdiss(kk), diss_corr.eta(kk), diss_corr.modStr(kk)]...
            = ModifiedStructureFun(stFun.Sdiss(kk).*normFac);

    end
    
    % Re-normalize the dissipation (this is a little silly isn't it...)
     diss_corr.St_diss = diss_corr.Rdiss ./ normFac;

    for jj = 1:length(diss_fields)
        % Take the mean of the gradient methods along each column of all jj
        % versions of the grad_diss
        diss_corr.(diss_fields{jj}).mGdiss = nanmean(diss_dat.(diss_fields{jj}));
        % calculate beta
        diss_corr.(diss_fields{jj}).beta = diss_corr.St_diss...
            ./ diss_corr.(diss_fields{jj}).mGdiss;
        % Modify the gradient values with beta 
        for kk = 1:length(diss_corr.(diss_fields{jj}).beta)

            diss_corr.(diss_fields{jj}).Mgrad(:,kk) = diss_corr.(diss_fields{jj}).beta(kk)...
                .* diss_dat.(diss_fields{jj})(:,kk);

        end

        % Take the mean of the modified gradient method at each column 
        diss_corr.(diss_fields{jj}).mMgrad = nanmean(diss_corr.(diss_fields{jj}).Mgrad);
        
    end
    
    diss_corr.Xmat = diss_dat.Xmat;
    diss_corr.Ymat = diss_dat.Ymat;
    
    % Saving stuff
    
    if ~exist(Sstr.strPth, 'dir')
        mkdir(Sstr.strPth)
    end
    
    if ~exist(Sstr.dissPth, 'dir')
        mkdir(Sstr.dissPth)
    end
    
    stSaveName = [Sstr.strPth filesep Lstr.names{ii} '_struc_funs'];
    dissSaveName = [Sstr.dissPth filesep Lstr.names{ii} '_diss_corr'];
    

    save(stSaveName, '-struct', 'stFun');
    save(dissSaveName, '-struct', 'diss_corr');
    
    clearvars -except Lstr Sstr problem fullTime sumTime Dj c
    
end

probSaveName = [Sstr.probPth filesep 'Problems'];
save(probSaveName, 'problem')

