tic
%% Parameter
%AllDataPath ='/home/michael/Cloud/PhD/MarsupialData/marsupial-data/ComparisonMaps/OnlineAnalysisOPM/TestData/';
AllDataPath = '~/Cloud/PhD/data/data share/Dunnart data/';


refWin=[1:10];
sigWin=[31:35];
stim_order = [0   22.5000   45.0000   67.5000   90.0000  112.5000  135.0000  157.5000 nan ];
lengthProngedElectrode_mm = 1.5;
pixels_per_mm = 20/0.48;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ADJUST FOR EACH ANIMAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ID = 'dunnartR';
expIds =  {[9 11 13 15],[10 12 14 16]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partId = ['A','B'];

%%% Powerprofile
profile_range_mm = [0.01 2];
profile_step_mm = 0.01;

%%% Pinwheel plateaus
lowpass_cutoffs = 0.15:0.01:1.5;

%%%Bootstrapping
Bootstrapsamples = 100;
alpha = 0.05;
apply_filter = true;

border = 10;

%% Data path
data_path = [AllDataPath ID '/'];
IntermedietResultsFilePW = [data_path 'PW_Stats.mat'];
IntermediatResultsFileOPM = [data_path 'OrientationStats.mat'];

%% Read in Data

BloodVesselImg = getBloodVesselImg(data_path,ID);
close

if ~isfile([data_path ID '.mat'])

    data = NoAveragePreProcessRawDataJason(expIds,refWin,sigWin,partId,data_path,ID);
    field_size_pix = size(data,1:2);

    ROI = getROI(data_path,ID);

    %% Make data info and data object

    data_info = makeDataInfoDunnart(ID,field_size_pix,pixels_per_mm,stim_order,expIds,refWin,sigWin,partId);
    data_obj = data_handle_corrected(data_info,data,ROI);

    %% Determine filter parameter from power spectrum

    power_profile = define_filter_settings(data_info,ROI,data,profile_range_mm,profile_step_mm);

    figure
    plot(power_profile.scale_mm,power_profile.values);
    xlabel('Scale in mm')
    ylabel('Power')
    xlim([0 max(profile_range_mm)])
    set(gca,'fontsize',15)

    prompt = "What is the lowpass filter value in mm? ";
    lowpass_mm = input(prompt);
    data_obj.set_filter_parameters('lowpass',lowpass_mm)
    data_info.settings.lowpass_mm = lowpass_mm;

    prompt = "What is the highpass filter value in mm? ";
    highpass_mm = input(prompt);
    data_obj.set_filter_parameters('highpass',highpass_mm)
    data_info.settings.highpass_mm = highpass_mm;

    close all


    %% determine lowpass filter from pinwheel plateau


    figure
    filtersPwNumber = find_lowpassPwNumber(data_obj,data_info,data_obj.filter_parameters.highpass,data_obj.filter_parameters.lowpass,lowpass_cutoffs);
    xlabel('lowpassfilter cutoff [mm]')
    ylabel('pinwheel number')

    prompt = "What is the lowpass filter value in mm based on the pinwheel density plateau? ";
    lowpass_mm = input(prompt);
    data_obj.set_filter_parameters('lowpass',lowpass_mm)
    data_info.settings.lowpass_mm = lowpass_mm;


    %% Last Check Filter parameter

    figure

    plot(power_profile.scale_mm,power_profile.values,'DisplayName','power profile unfiltered map');
    hold on
    plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','lowpass cutoff')
    hold on
    plot([data_obj.filter_parameters.highpass data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','highpass cutoff')

    xlabel('Scale in mm')
    ylabel('Power')
    xlim([0.01 2])
    ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
    set(gca,'fontsize',15)
    legend()


    disp(['lowpassfilter cutoff [mm]: ' num2str(data_obj.filter_parameters.lowpass)])
    disp(['highpassfilter cutoff [mm]: ' num2str(data_obj.filter_parameters.highpass)])

    prompt = "Are the filter parameter right?  ";
    Answer = char(lower(input(prompt)));

    if Answer(1) ~= 'y' 
        prompt = "What is the lowpass filter value in mm? ";
        lowpass_mm = input(prompt);
        data_obj.set_filter_parameters('lowpass',lowpass_mm)
        data_info.settings.lowpass_mm = lowpass_mm;

        prompt = "What is the highpass filter value in mm? ";
        highpass_mm = input(prompt);
        data_obj.set_filter_parameters('highpass',highpass_mm)
        data_info.settings.highpass_mm = highpass_mm;
    end

    close all

    %% save (filter)parameter and data
    save_info(data_path,data_info)
    save([data_path ID],'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
else
    %% load filter data
    load([data_path ID '.mat'],'data_obj','data_info','lowpass_cutoffs','filtersPwNumber','power_profile')
    
    %% plot pinwheel plateau
    
    figure()
    plot(filtersPwNumber.global_plateau.lowpass_vs_density(:,1),filtersPwNumber.global_plateau.lowpass_vs_density(:,2))
    xlabel('lowpassfilter cutoff [mm]')
    ylabel('pinwheel number')
    
    
    %% plot filter plus power spectrum
    figure()

    plot(power_profile.scale_mm,power_profile.values,'DisplayName','power profile unfiltered map');
    hold on
    plot([data_obj.filter_parameters.lowpass data_obj.filter_parameters.lowpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','lowpass cutoff')
    hold on
    plot([data_obj.filter_parameters.highpass data_obj.filter_parameters.highpass],[min(power_profile.values,[],'all') max(power_profile.values,[],'all')],'DisplayName','highpass cutoff')

    xlabel('Scale in mm')
    ylabel('Power')
    xlim([0.01 2])
    ylim([min(power_profile.values,[],'all') max(power_profile.values,[],'all')])
    set(gca,'fontsize',15)
    legend()


    disp(['lowpassfilter cutoff [mm]: ' num2str(data_obj.filter_parameters.lowpass)])
    disp(['highpassfilter cutoff [mm]: ' num2str(data_obj.filter_parameters.highpass)])
end

%% Prepare Certainty analysis

data_obj.prepare_samples_array(Bootstrapsamples)

%% get PW certainty and confidence intervals



if isfile(IntermedietResultsFilePW)
    disp('PW Stats already exist')
    load(IntermedietResultsFilePW,'pinwheel_stats','data_obj','data_info')
else
    tracker_obj = pinwheel_tracker;
    simple_track=true;
    [pinwheel_stats,pinwheel_spurious] = get_pinwheel_stats(data_obj,tracker_obj,simple_track);
    save(IntermedietResultsFilePW,'pinwheel_stats','pinwheel_spurious')
end
    


% get OPM selectivity and confidence


if isfile(IntermediatResultsFileOPM)
    disp('Orientation Stats already exist')
    load(IntermediatResultsFileOPM,'orientation_stats');
else
    orientation_stats = get_orientation_stats(data_obj,alpha,apply_filter);
    save(IntermediatResultsFileOPM,'orientation_stats','data_obj')
end


%% Plot results
ElectrodePositions0 = {[border,border,border+lengthProngedElectrode_mm*pixels_per_mm,border]};
PlotMapsForElectrodePositions(pinwheel_stats,orientation_stats,ElectrodePositions0,data_obj,data_info,BloodVesselImg,data_path)

%% Input Electrode positions
N_ElectrodePronges = 1;
NElectrodes = 0;
ElectrodePositions = {};
while N_ElectrodePronges ~= 0
    
    prompt = "How many Prongs should the next electrode have? Options are 1 or 4!  ";
    N_ElectrodePronges = input(prompt);
    
    if N_ElectrodePronges == 1
        prompt = 'Input the Electrode pixel position in the format [x y]!  ';
        answer = input(prompt);
        if size(answer,2) == 2
            NElectrodes = NElectrodes+1;
            ElectrodePositions{NElectrodes}= answer;
        end
    elseif N_ElectrodePronges == 4
        prompt = 'Input the Electrode pixel position in the format [x1 y1 x2 y2]!  ';
        answer = input(prompt);
        if size(answer,2) == 4
            NElectrodes = NElectrodes+1;
            ElectrodePositions{NElectrodes}= AdjustElectrodeWidth(answer,lengthProngedElectrode_mm*pixels_per_mm);
        end
    else
        if N_ElectrodePronges ~=0
            disp(['An electrode with ' num2str(N_ElectrodePronges) ' pronges does not exist!  '])
        end
    end
    
end

close all

%% Plot Maps with defined Electrodes
PlotMapsForElectrodePositions(pinwheel_stats,orientation_stats,ElectrodePositions,data_obj,data_info,BloodVesselImg,data_path)


%% Adjust Electrodes
AdjustElectrode = 1;

while AdjustElectrode ~= 0
    prompt = ['Should any of the ' num2str(size(ElectrodePositions,2)) ' electrode be adjusted? If yes, enter its number! If not enter 0!  '];
    AdjustElectrode = input(prompt);
    if AdjustElectrode ~=0
        if AdjustElectrode <= size(ElectrodePositions,2)
            prompt = ['Should the size of electrode ' num2str(AdjustElectrode) ' be adjusted? If yes, enter the new size! If not enter 0! Enter negative number to delete Electrode!  '];
            Answer = input(prompt);
            if Answer == 0
                N_ElectrodePronges = (size(ElectrodePositions{AdjustElectrode},2)-2)*(3/2)+1;
            else
                N_ElectrodePronges = Answer;
            end
            disp('Previous position:')
            disp(ElectrodePositions{AdjustElectrode})
        else
            AdjustElectrode = size(ElectrodePositions,2)+1;
            prompt = ['You want to add a new electrode? Please enter its number of pronges! If not answer 0!  '];
            N_ElectrodePronges = input(prompt);
        end
        if N_ElectrodePronges == 1
            prompt = 'Input the Electrode pixel position in the format [x y]!  ';
            answer = input(prompt);
            if size(answer,2) == 2
                ElectrodePositions{AdjustElectrode}= answer;
            end
        elseif N_ElectrodePronges == 4
            prompt = 'Input the Electrode pixel position in the format [x1 y1 x2 y2]!  ';
            answer = input(prompt);
            if size(answer,2) == 4
                ElectrodePositions{AdjustElectrode}= AdjustElectrodeWidth(answer,lengthProngedElectrode_mm*pixels_per_mm);
            end
        elseif N_ElectrodePronges < 0
            ElectrodePositions{AdjustElectrode}= [];
            disp(['Former electrode ' num2str(AdjustElectrode) ' has been deleted.'])
        else
            if N_ElectrodePronges ~=0
                disp(['An electrode with ' num2str(N_ElectrodePronges) ' pronges does not exist!  '])
            end
        end
    end
    close all
    PlotMapsForElectrodePositions(pinwheel_stats,orientation_stats,ElectrodePositions,data_obj,data_info,BloodVesselImg,data_path)
end

%% Plot and save final Blood vessel map

figure()
PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)),1)
PlotColoredElectrodes(ElectrodePositions)
savefig([data_path 'BloodVesselsWithElectrodes'])

%% Disp elapsed time
disp(toc/60)


%% Functions
function data_info = makeDataInfoDunnart(ID,field_size_pix,pixels_per_mm,stim_order,expIds,refWin,sigWin,partId)

    clear data_info

    data_info.ID = ID;
    data_info.animal = 'Dunnart';
    data_info.field_size_pix = field_size_pix;
    data_info.pix_per_mm = pixels_per_mm;
    data_info.pixels_per_mm = data_info.pix_per_mm;
    data_info.stim_order = stim_order;
    data_info.expIds = expIds;
    data_info.refWin = refWin;
    data_info.sigWin = sigWin;
    data_info.partId = partId;

end

function BloodVesselImg = getBloodVesselImg(data_path,ID)
    fig = openfig([data_path ID '.fig']);
    BloodVesselImg = getimage(fig);
end

function ROI = getROI(data_path,ID)
    load([data_path ID '_Mask.mat'],'a')
    ROI = a;
end
