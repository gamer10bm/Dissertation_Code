function [out_struct,failcheck_mat,abs_min_off_m,abs_min_t] = constellation_calcs(...
    sats,driving_sat_id,avg_array_rangem,fail_min_off_m,fail_prob)

Settings_Targets
%Parse inputs
num_states = length(sats(driving_sat_id).states);
num_sats = length(sats);
avg_array_maxm = max(avg_array_rangem);
avg_array_minm = min(avg_array_rangem);
non_drive_sats = 1:num_sats;
non_drive_sats(driving_sat_id) = [];

%% Do post calculations for constellation relative positions
%Initialize
state_pos_cells = cell(1,num_states);
alt_mat = zeros(num_states,num_sats);
along_pos_mat = zeros(num_states,num_sats);
cross_pos_mat = zeros(num_states,num_sats);
elev_pos_mat = zeros(num_states,num_sats);
inplanemag_pos_mat = zeros(num_states,num_sats);
mag_pos_mat = zeros(num_states,num_sats);
avg_mag_m = zeros(1,num_states);
avg_inplanemag_m = zeros(1,num_states);
Grecheck = zeros(1,num_states);
Antcheck = zeros(1,num_states);
datacheck = zeros(1,num_states);
datametric_totallinekm = zeros(1,num_states); %line-km of collected data total
datametric_Grelinekm = zeros(1,num_states); %line-km of data from Greenland
datametric_Antlinekm = zeros(1,num_states); %line-km of data from Antarctica

failcheck_mat = zeros(num_states,num_sats);
abs_min_off_m = 1e8;
abs_min_t = 0;
%Iterate
for sid = 1:num_states
    %% Convert all positions to RSW and calculate
    %Initialize short calcs
    sat_pos_cells = cell(1,num_sats);
    alt_vec = zeros(1,num_sats);
    along_pos = zeros(1,num_sats);
    cross_pos = zeros(1,num_sats);
    elev_pos = zeros(1,num_sats);
    inplanemag_pos = zeros(1,num_sats);
    mag_pos = zeros(1,num_sats);
    %Rotate current positions
    for sat_id = 1:num_sats
        %Grab position vectors
        satnow_Rvec = sats(sat_id).states(sid).R;
        %Grab orbital elements
        satnow_OE = sats(sat_id).states(sid).OE;
        Omega = satnow_OE.RA_RightAscension;
        argper = satnow_OE.w_ArgumentofPerigee;
        incl = satnow_OE.i_Inclination;
        %Calculate rotation matrix
        [Rmat] = ECI2ConstFrame(satnow_OE);
        %Rotate position vector
        sat_pos_cells{sat_id} = Rmat*satnow_Rvec;        
    end
    %Calculate relative positions based on driving sat
    satdrive_pos = sat_pos_cells{driving_sat_id};
    for sat_id = 1:num_sats
        satnow_pos = sat_pos_cells{sat_id};
        alt_vec(sat_id) = sats(sat_id).states(sid).alt;
        along_pos(sat_id) = satdrive_pos(1) - satnow_pos(1);
        cross_pos(sat_id) = satdrive_pos(2) - satnow_pos(2);
        elev_pos(sat_id) = satdrive_pos(3) - satnow_pos(3);
        inplanemag_pos(sat_id) = sqrt((cross_pos(sat_id))^2+(elev_pos(sat_id))^2);
        mag_pos(sat_id) = sqrt((along_pos(sat_id))^2+(cross_pos(sat_id))^2+(elev_pos(sat_id))^2);
    end
    
    non_drive_satids = [1:driving_sat_id-1,driving_sat_id+1:num_sats];
    %Calculate averages
    avg_inplanemag_m(sid) = mean(abs(diff(inplanemag_pos(non_drive_satids))))*1e3;
    avg_mag_m(sid) = mean(abs(diff(mag_pos(non_drive_satids))))*1e3;
    
    %Do failure check
    for sat_id = 1:num_sats
        satcheck_pos = sat_pos_cells{sat_id};
        %Check for previous failure
        if sid == 1 || (sat_id<num_sats && ~failcheck_mat(sid-1,sat_id))
            for oid = sat_id+1:num_sats
                satother_pos = sat_pos_cells{oid};
                check_mag_m = norm(satcheck_pos-satother_pos)*1000;
                %Update abs_min_off_m
                if check_mag_m<abs_min_off_m
                    abs_min_off_m = check_mag_m;
                    abs_min_t = sats(sat_id).states(sid).t;
                end
                %Determine if failure
                if check_mag_m < fail_min_off_m
                    %Roll for failure save
                    randroll = rand();
                    if randroll <= fail_prob
                        %Failure has occured between these two sats
                        failcheck_mat(sid,sat_id) = 1;
                        failcheck_mat(sid,oid) = 1;
                    else
                        %Probability saved the day
                        fprintf('\tThat was close! %s almost hit %s\n',sats(sat_id).name,sats(oid).name)
                    end
                end
            end
        else
            %Mark previous failure
            failcheck_mat(sid,sat_id) = failcheck_mat(sid-1,sat_id);
        end
    end
    
    %Do checks for radar targets
    satdrive_lat = sats(sat_id).states(sid).lat;
    satdrive_lon = sats(sat_id).states(sid).lon;
    Grecheck(sid) = abs(geodist(radar_Gre.Lat,radar_Gre.Lon,satdrive_lat,satdrive_lon))<radar_Gre.Rad;
    if ~Grecheck(sid)
        Antcheck(sid) = abs(geodist(radar_Ant.Lat,radar_Ant.Lon,satdrive_lat,satdrive_lon))<radar_Ant.Rad;
    else
        Antcheck(sid) = false;
    end
    
    %Collect data
    if sid>1 && (Grecheck(sid) || Antcheck(sid)) && (avg_inplanemag_m(sid)<avg_array_maxm && avg_inplanemag_m(sid)>avg_array_minm)
        datacheck(sid) = true;
        %Check if segment of data collected
        if datacheck(sid-1)
            %Accumulate line-km of data
            segdist_km = abs(geodist(satdrive_lat,satdrive_lon,...
                sats(sat_id).states(sid-1).lat,sats(sat_id).states(sid-1).lon));
            datametric_totallinekm(sid) = datametric_totallinekm(sid-1) + segdist_km;
            if Grecheck(sid)
                datametric_Grelinekm(sid) = datametric_Grelinekm(sid-1) + segdist_km;
            elseif Antcheck(sid)
                datametric_Antlinekm(sid) = datametric_Antlinekm(sid-1) + segdist_km;
            end
        end
    end
    %Fix metric accumulation
    if datametric_totallinekm(sid) == 0 && sid>1
        datametric_totallinekm(sid) = datametric_totallinekm(sid-1);
    end
    if datametric_Grelinekm(sid) == 0 && sid>1
        datametric_Grelinekm(sid) = datametric_Grelinekm(sid-1);
    end
    if datametric_Antlinekm(sid) == 0 && sid>1
        datametric_Antlinekm(sid) = datametric_Antlinekm(sid-1);
    end
    
    
    %Load to outputs
    state_pos_cells{sid} = sat_pos_cells;
    alt_mat(sid,:) = alt_vec;
    along_pos_mat(sid,:) = along_pos;
    cross_pos_mat(sid,:) = cross_pos;
    elev_pos_mat(sid,:) = elev_pos;
    inplanemag_pos_mat(sid,:) = inplanemag_pos;
    mag_pos_mat(sid,:) = mag_pos;
end

%Load output structure
out_struct = struct('state_pos',{state_pos_cells},'alt',alt_mat,...
    'along_pos',along_pos_mat,'cross_pos',cross_pos_mat,'elev_pos',elev_pos_mat,...
    'inplane_mag',inplanemag_pos_mat,'total_mag',mag_pos_mat,'avg_total_mag',avg_mag_m,...
    'avg_inplane_mag',avg_inplanemag_m,'Greencheck',Grecheck,'Antcheck',Antcheck,...
    'datacheck',datacheck,'data_total',datametric_totallinekm,...
    'data_Green',datametric_Grelinekm,'data_Ant',datametric_Antlinekm);