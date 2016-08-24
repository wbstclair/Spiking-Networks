function data = networkEvaluate(currentIDDM, currentIDDV, currentEDDM, currentEDDV, currentINCPN, currentENCPN, currentStimMean, currentStimVariance, rateCode, stimulationKind, timeThreshold, returnSuperData)
% delay distribution (internally) within a cluster. mean and variance
%internalDelayDistribMean = [1:6:30]; IDDM
%internalDelayDistribVariance = [0:4:18]; IDDV
% delay between clusters (externally) mean and variance
%externalDelayDistribMean = [1:6:30]; EDDM
%externalDelayDistribVariance = [0:4:18]; EDDV
% connections per neuron within a cluster (corresponds to density)
%internalNumConnectionsPerNeuron = [10:20:90]; INCPN
% connections per neuron between clusters (corresponds to density)
%externalNumConnectionsPerNeuron = [10:10:50]; ENCPN
% stimuli firing rate (per second) mean and variance
%stimuliMean = [30:100:500]; StimMean
%stimuliVariance = [0:10:40]; StimVariance % irrelevant but needed if
%poisson (for arguments to match)
%rateCode = rand(200,1); row # = neuron index, value = rate
%stimulationKind = 1, 2. Stimulation kind 1 implies poisson. Stimulation
%kind 2 implies the (fixed) precise frequency+jitter method.
%timeThreshold , used to impose a hard maximum on how long it takes to run
%the evaluation. The number itself refers to the # of spikes in the list,
%not the clock time it takes to run the evaluation.
% Return superData, which has membrane potential data and network
% connectivity.
% VERSION 5.0, or so named.
%
% This code was written by William Benjamin St. Clair from 2013-2016
% to generate the parametric analysis used in "Contexts for the Emergence of
% Significant Neural Dynamics": http://escholarship.org/uc/item/28v2w6rz 


if ~exist('timeThreshold')
    timeThreshold = NaN;
end
if ~exist('returnSuperData')
    returnSuperData = 0;
end

rng('shuffle');

numClusters = 6;


numNuclei = numClusters;
% Initialize the metanetwork cell array.
%nuclei = cell(numNuclei,numNuclei);
for i = 1:numNuclei
    for j = 1:numNuclei
        nuclei{i}{j} = [];
    end
end


% network parameters.
%Number of excitatory neurons
params.Ne = 200;
%Number of inhibitory neurons
params.Ni = 40;
params.N = params.Ne + params.Ni;
% neuron model parameters for the excitatory neurons
params.exA = 0.02;
params.exB = 0.2;
params.exC = -65;
params.exD = 8;
% neuron model parameters for the inhibitory neurons.
params.inA = 0.1;
params.inB = 0.2;
params.inC = -65;
params.inD = 2;
% Number of outgoing connections per neuron. 
% I have set this to zero to allow no self-connections for the source
% nuclei.
params.numConnectionsPerNeuron = 0;
% a vector of inclusive lower bound and inclusive upper bound of the range of conductance delays. in ms
% Izhi&Szat have 20ms max, uniform.
params.delayRange = currentIDDM + [0 currentIDDV];
%params.isSpatialNetwork = true;
%params.spatialClusterType = 'poisson';
params.initExWeight = 5;
params.initInWeight = -4;
% weight bounds, applied every millisecond.
params.weightUpperbound = 8;
params.weightLowerbound = -8;
% Izhikevich iters two .5 ms to add to 1ms. This is not really a parameter.
params.timeStep = .5; % 0.5 ms.

network1 = networkBuild(params);
network1.nucleiIndex = 1;

% network parameters.
%Number of excitatory neurons
params.Ne = 200;
%Number of inhibitory neurons
params.Ni = 40;
params.N = params.Ne + params.Ni;
% neuron model parameters for the excitatory neurons
params.exA = 0.02;
params.exB = 0.2;
params.exC = -65;
params.exD = 8;
% neuron model parameters for the inhibitory neurons.
params.inA = 0.1;
params.inB = 0.2;
params.inC = -65;
params.inD = 2;
% Number of outgoing connections per neuron.
params.numConnectionsPerNeuron = currentINCPN;
% a vector of inclusive lower bound and inclusive upper bound of the range of conductance delays. in ms
% Izhi&Szat have 20ms max, uniform.
params.delayRange = currentIDDM + [0 currentIDDV];
params.initExWeight = 5;
params.initInWeight = -4;
% weight bounds, applied every millisecond.
params.weightUpperbound = 8;
params.weightLowerbound = -8;
% Izhikevich iters two .5 ms to add to 1ms. This is not really a parameter.
params.timeStep = .5; % 0.5 ms.


nuclei{1}{1} = network1;
clear network1;

% Build numClusters # of nuclei.
for currentCluster = 2:numClusters
    network = networkBuild(params);
    network.nucleiIndex = currentCluster;
    nuclei{currentCluster}{currentCluster} = network;
end
clear params network currentCluster


% Now, build the connections between the source nuclei and the different
% clusters.
delayRange = currentEDDM + [0 currentEDDV];
%baseDelay = 5; % the minimum delay
%delayVariability = 3;
numConnectionsPerNeuron = currentENCPN;
weightUpperbound = 8;
weightLowerbound = -8;
numNuclei = numClusters;

% this builds connections between the source cluster and downstream
% clusters (forward only)
for currentCluster = 2:numClusters
    % Feed forward connect nuclei index to another nuclei index.
    ffConnect = [1 currentCluster];
    % a vector of inclusive lower bound and inclusive upper bound of the range of conductance delays. in ms
    nuclei = buildConnections(nuclei, numNuclei, ffConnect, delayRange, ...
        numConnectionsPerNeuron, weightLowerbound, weightUpperbound);
end

% This builds connections between each downstream cluster.
for sourceCluster = 2:numClusters
    for targetCluster = 2:numClusters
        % Only build external connections.
        if sourceCluster ~= targetCluster
            ffConnect = [sourceCluster targetCluster];
            nuclei = buildConnections(nuclei, numNuclei, ffConnect, delayRange, ...
              numConnectionsPerNeuron, weightLowerbound, weightUpperbound);
        end
    end
end
clear ffConnect numConnectionsPerNeuron
clear weightUpperbound weightLowerbound delayRange        

%Init stimuli. It should yield a cell stim with stimNumber of elements, each
%referring to a 1000ms period. Each spike stimuli is a row element of
%stim{stimNumber}, consisting of [neuron timeOfSpike stimmedNuclei]


% a vector pattern of stimulation, where 0 implies no stimulus.
stimPattern = [1];
% setting runForWholeStimPattern to anything besides 1 will only allow
%   the network to run for the first stimuli.
runForWholeStimPattern = 1;

stimLength = 1000;

% Use the provided rate code, assume neuron index implied by the rate code index. 
stimFrequencies = currentStimMean .* rateCode(:,1);
%stimFrequencies becomes the Hz for that particular neuron.
stimISIs = stimLength./stimFrequencies;

switch stimulationKind
    case 1 % poisson distribution
        fullStimuli = [];
        for i = 1:length(stimFrequencies)
            if stimISIs(i) > 0 
                neuronTimes = poissrnd(stimISIs(i));
                while neuronTimes(end)-neuronTimes(1) < stimLength
                    neuronTimes = [neuronTimes; neuronTimes(end) + poissrnd(stimISIs(i))];
                end
                % In order for the stimulus to be properly centered, we create the
                % spikes until 1000ms has elapsed from the very first spike in the
                % code. Since we don't know which neuron will be the first neuron
                % to spike, we generate all the spikes first, then prune and shift
                % once everything has been computed.
                neuronTimes = neuronTimes(neuronTimes < (stimLength+neuronTimes(1)));     
                fullStimuli = [fullStimuli; neuronTimes i*ones(length(neuronTimes),1)];
            end
        end
        fullStimuli(:,1) = fullStimuli(:,1) - min(fullStimuli(:,1)) + 1;
        fullStimuli = fullStimuli(fullStimuli(:,1) < 1000,:);
        stim{1} = [fullStimuli ones(size(fullStimuli,1),1)];
    case 2 % normal distribution
        fullStimuli = [];
        for i = 1:length(stimFrequencies)
            if stimISIs(i) > 0 
                neuronTimes = round(normrnd(stimISIs(i),currentStimVariance));
                if neuronTimes <= 0
                    neuronTimes = 1;
                end
                while neuronTimes(end)-neuronTimes(1) < stimLength
                    newTime = round(normrnd(stimISIs(i),currentStimVariance));
                    if newTime <= 0
                        newTime = 1;
                    end
                    neuronTimes = [neuronTimes; neuronTimes(end) + newTime];
                end
                % In order for the stimulus to be properly centered, we create the
                % spikes until 1000ms has elapsed from the very first spike in the
                % code. Since we don't know which neuron will be the first neuron
                % to spike, we generate all the spikes first, then prune and shift
                % once everything has been computed.
                neuronTimes = neuronTimes(neuronTimes < (stimLength+neuronTimes(1)));
                neuronTimes = round(neuronTimes);
                fullStimuli = [fullStimuli; neuronTimes i*ones(length(neuronTimes),1)];
            end
        end
        fullStimuli(:,1) = fullStimuli(:,1) - min(fullStimuli(:,1)) + 1;
        fullStimuli = fullStimuli(fullStimuli(:,1) < 1000,:);
        stim{1} = [fullStimuli ones(size(fullStimuli,1),1)];
end



% This was how stimuli were created with jitter independently varied.
% % First column is the neuron number, the second column is its spiking in Hz
% stimFrequencies = ceil(currentStimMean * rateCode(:,1));
% 
% fullStimuli = [];
% for i = 1:length(rateCode)
% 
%     % How often, in milliseconds, the neuron in question is to be
%     % stimulated.
%     stimFreq = floor(1000/stimFrequencies(i));
% 
% 	for preciseStimNumber = 1:stimFrequencies(i)
%         preciselyVariedTime = round(normrnd(1000/stimFrequencies(i)*preciseStimNumber, currentStimVariance));
%         if preciselyVariedTime <= stimLength && preciselyVariedTime >= 1
%             fullStimuli = [fullStimuli; preciselyVariedTime i 1];
%         end
%     end
% end
% stim{1} = fullStimuli;
% clear fullStimuli stimFreq


% Determine background noise. Must be between 0 and 1.
% This corresponds to the likelihood that a neuron will fire due to noise
% per second.
bgNoise = 0; %.05;
% The amplitude of the noise, in mV.
noiseAmplitude = 20;

% Since there will not be any continued runs here, totalTime = 0 (always).
totalTime = 0;
firings = [];
%'running the network...'
previousTime = totalTime;
if runForWholeStimPattern == 1
    runTime = length(stimPattern)*stimLength;
else
    runTime = stimLength;
end
if returnSuperData
    superData = cell(runTime,1);
end
% Run the network for the entirety of runTime.
%fprintf('\n');
for t = (previousTime+1):(previousTime+runTime)
    totalTime = t;
    if stimPattern(ceil((t - previousTime)/1000)) == 0
        stimulusPattern = [];
    else
        stimulusPattern = stim{stimPattern(ceil((t - previousTime)/1000))};
    end
    for nucleiNumber = 1:numNuclei
%    for nucleiNumber = numNuclei:-1:1
        % Once per second, for each nuclei, calculate the background noise.
		% TODO: Change to only a ms!
        willFire = find(rand(nuclei{nucleiNumber}{nucleiNumber}.params.N,1) < bgNoise/1000);
        %neuronFireTimes = ceil(rand(length(willFire),1)*1000) + previousTime;
        %inputList = [neuronFireTimes willFire];
%        x= [0 0 0 0 0];
        networkIterate
%        allTocs(t,:) = x;
    end
%     if mod(t, 100) == 0
%         fprintf([num2str(t) ' ']);
%         fprintf([num2str(length(firings)) ' ']);
%     end
%     if length(firings) > timeThreshold
%         break;
%     end
    if returnSuperData
        superData{t} = nuclei;
    end

end
clear willFire neuronFireTimes nucleiNumber previousTime t

%data = firings;
if returnSuperData
    data.firings = firings;
    data.superData = superData;
else
    data = firings;
end

return