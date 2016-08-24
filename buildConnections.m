function nuclei = buildConnections(nuclei, numNuclei, ffConnect, ...
    delayRange, numConnectionsPerNeuron, weightLowerbound, weightUpperbound, allowInhibProjections)
% function wholeNetwork = buildConnections(nuclei, numNuclei, ffConnect, ...
%     baseDelay, numConnectionsPerNeuron, weightLowerbound, weightUpperbound)
%
% build feed forward connections from nuclei ffConnect(1) to ffConnect(2)
% with the baseDelay minimum conductance delay and numConnectionsPerNeuron
% connections per neuron in nuclei ffConnect(1).
%
% weightLowerbound is the lower bound to be enforced for the connection weights,
% 		   corresponds to the inhibitory minimum.
% weightUpperbound is the upper bound to be enforced for the connection weights,
% 		   corresponds to the excitatory maximum.

% connectivityMatrix is a matrix of zeros and ones, with ones implying the
% existence of a connection. It is used to impose the network connectivity
% on any matrix of equal size.
%
% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012

if ~exist('allowInhibProjections')
    allowInhibProjections = 0;
end

params1 = nuclei{ffConnect(1)}{ffConnect(1)}.params;
params2 = nuclei{ffConnect(2)}{ffConnect(2)}.params;

% Only make inhibitory projections 
if allowInhibProjections
    connectivityMatrix = zeros(params1.N, params2.N);
    for i = 1:params1.N
        connectionsSoFar = [];
        for j = 1:numConnectionsPerNeuron
            if i <= params1.Ne
                randomConnection = ceil(rand*params2.N);
            else
                randomConnection = ceil(rand*params2.Ne);
            end
            while nnz(connectionsSoFar == randomConnection) > 0
                if i <= params1.Ne
                    randomConnection = ceil(rand*params2.N);
                else
                    randomConnection = ceil(rand*params2.Ne);
                end
            end
            connectionsSoFar = [connectionsSoFar randomConnection];
            connectivityMatrix(i,randomConnection) = 1;
        end
    end
    clear randomConnection
    clear connectionsSoFar
else
    connectivityMatrix = zeros(params1.N, params2.N);
    for i = 1:params1.Ne
        connectionsSoFar = [];
        for j = 1:numConnectionsPerNeuron
            if i <= params1.Ne
                randomConnection = ceil(rand*params2.N);
            else
                randomConnection = ceil(rand*params2.Ne);
            end
            while nnz(connectionsSoFar == randomConnection) > 0
                if i <= params1.Ne
                    randomConnection = ceil(rand*params2.N);
                else
                    randomConnection = ceil(rand*params2.Ne);
                end
            end
            connectionsSoFar = [connectionsSoFar randomConnection];
            connectivityMatrix(i,randomConnection) = 1;
        end
    end
    clear randomConnection
    clear connectionsSoFar
end

% with uniform weights.
S=[params1.initExWeight*ones(params1.Ne,params2.N); ...
   params1.initInWeight*ones(params1.Ni,params2.N)];
S = S .* connectivityMatrix;

%Izhi-Szat paper has short-term STDP values which scale input for each
%neuron, unique to each synapse.
%    gatingValues = zeros(size(S));

% The millisecond delays between spikes.
% The range of the millisecond delays between spikes, uniformly chosen.
delayRange = delayRange(1):delayRange(2);
conductanceDelays = delayRange(ceil(rand(params1.N,params2.N)*length(delayRange)));

% this line normally sets all inhibitory connections to delay 1, as in the
% izhikevich paper. However, for nuclei to nuclei we expect atleast a
% minimum delay of baseDelay+1+1 (the original 1, but for each nuclei).
if params1.N > params1.Ne
	conductanceDelays((params1.Ne+1):params1.N,:) = delayRange(1);
end
conductanceDelays = conductanceDelays .* connectivityMatrix;







% Initially, there are NO conducting potentials (assumption).
conductingPotentials = [];

nuclei{ffConnect(1)}{ffConnect(2)}.conductingPotentials = conductingPotentials;
nuclei{ffConnect(1)}{ffConnect(2)}.conductanceDelays = conductanceDelays;
nuclei{ffConnect(1)}{ffConnect(2)}.S = S;
nuclei{ffConnect(1)}{ffConnect(2)}.connectivityMatrix = connectivityMatrix;
nuclei{ffConnect(1)}{ffConnect(2)}.lastFire = ...
    NaN(nuclei{ffConnect(1)}{ffConnect(1)}.params.N, ...
    nuclei{ffConnect(2)}{ffConnect(2)}.params.N);
nuclei{ffConnect(1)}{ffConnect(2)}.params.weightLowerbound = weightLowerbound;
nuclei{ffConnect(1)}{ffConnect(2)}.params.weightUpperbound = weightUpperbound;

% init a lastFire = [] for each empty region.

for n1 = 1:numNuclei
	for n2 = 1:numNuclei
		if n1 ~= n2 & isempty(nuclei{n1}{n2})
			nuclei{n1}{n2}.lastFire = [];
		end
	end
end


end
