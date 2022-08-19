function [] = optimize_func(const_vec,OE0,nuangs,iangs,RAangs,wangs,init_utcvec,...
    Norbits,timesteps,avg_array_range_m,fail_dist_m,fail_prob,savedatapath,fnfrmt,wid,nuid,iid,RAid)


%% Initialize orbit
%Grab driving orbital elements
hmag = OE0(1); emag = OE0(2); iang = OE0(3); nuang = OE0(6); RAang = OE0(4); wang = OE0(5)-5*pi/100;
%Iterate over defined constellation
RXnum = 0; TXnum = 0;
for satid = 1:length(const_vec)
    %% Satellite Details
    %Grab constellation details
    satnow_str = const_vec{satid};
    %Name satellite
    if strcmp(satnow_str,'RX')
        satname = sprintf('%s%d',satnow_str,RXnum);
        RXnum = RXnum+1;
        %Add RX satellite
    elseif strcmp(satnow_str,'TX')
        satname = sprintf('%s%d',satnow_str,TXnum);
        TXnum = TXnum+1;
        %Add TX satellite
    else
        satname = satnow_str;
        %Add rando sat
    end

    %Initialize satellite model
    SC_model_fh = str2func(sprintf('initiate_%s_model',satnow_str));
    SCnow = SC_model_fh();

    %Take direct changes to orbital elements
    nuoff = nuangs(satid); ioff = iangs(satid); RAoff = RAangs(satid); woff = wangs(satid);

    %Calculate position and velocity vectors    
    [~,Rnow,Vnow] = coe2RV(hmag,emag,iang+ioff,nuang+nuoff,RAang+RAoff,wang+woff,SCnow.mu);
    COEstruct0 = RV2coe(Rnow,Vnow,SCnow.mu);

    % Load sats structure
    sats(satid).name = satname;
    sats(satid).SC = SCnow;
    sats(satid).Rinit = Rnow;
    sats(satid).Vinit = Vnow;
    sats(satid).COE0 = COEstruct0;
end
%% Propagate Orbit
[sats] = Sim_Constellation_matlab_J2Drag(sats,init_utcvec,Norbits,timesteps);

%% Do Post Calculations
driving_sat_id = find(strcmp('TX',const_vec)==1,1,'first');
num_states = length(sats(1).states);

[out_struct,failcheck_mat,minnow,mintnow] = constellation_calcs(sats,driving_sat_id,avg_array_range_m,fail_dist_m,fail_prob);

%% Save relavant data
%Send to save data
savefn = feval(fnfrmt,wid,nuid,iid,RAid);
save(fullfile(savedatapath,savefn),'out_struct','failcheck_mat','minnow','mintnow','nuangs','iangs','RAangs','wangs');
% fprintf('\n\t\t%s Saved\n',savefn);