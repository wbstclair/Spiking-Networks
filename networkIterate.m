% assume network is a cell 'nuclei'
% with a number of nuclei equally 'numNuclei'
%
%membranePotential =>
%       wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential
%
% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012





% Apply input if the stimulus pattern is nonempty for this nuclei.
if ~isempty(stimulusPattern(:,3) == nucleiNumber)
    stimulus = stimulusPattern(:,1) == mod(t-previousTime,1000) & ...
        stimulusPattern(:,3) == nucleiNumber;
    nuclei{nucleiNumber}{nucleiNumber}.membranePotential(stimulusPattern(stimulus,2)) = 30;
end

%%%% Find action potentials.
fired=find(nuclei{nucleiNumber}{nucleiNumber}.membranePotential>=30);    % indices of spikes

% For observation, we record the occurrence of the firings.
firings=[firings; t+0*fired,fired, nucleiNumber+0*fired];

% Instantiate threshold effects, setting membrane potential to c and
% adding d to the equalibrium force.
nuclei{nucleiNumber}{nucleiNumber}.membranePotential(fired)=nuclei{nucleiNumber}{nucleiNumber}.params.c(fired);
nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce(fired)=nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce(fired)+nuclei{nucleiNumber}{nucleiNumber}.params.d(fired);

% Since the Polychronization paper uses nonplastic inhibitory
% connections, we only should change weight values for excitatory
% neurons that fired. Hence, we obtain the excitatory subset of the
% fired neurons for this purpose.
excitatorySubset = fired(fired <= nuclei{nucleiNumber}{nucleiNumber}.params.Ne);
%x(1) = toc;
%tic
% Calculate STDP for connections in the nuclei based on new firings.
% Calculate STDP for all outgoing connections.
for receivingNuclei = 1:numNuclei
	% Only attempt to calculate STDP when there exists some connectivity!
	if ~isempty(nuclei{nucleiNumber}{receivingNuclei}.lastFire) ...
	 || ~isempty(nuclei{receivingNuclei}{nucleiNumber}.lastFire)
		% STDP applies spike time dependent plasticity
        % Uncomment this next line to disable STDP
		[stdpEffectsPos stdpEffectsNeg] = ...
			STDP(nuclei{nucleiNumber}{receivingNuclei}.lastFire, ...
			excitatorySubset, ...
			t, ...
			nuclei{receivingNuclei}{receivingNuclei}.params.Ne, ...
			nuclei{receivingNuclei}{nucleiNumber}.lastFire);

        if ~isempty(nuclei{receivingNuclei}{nucleiNumber}.lastFire)
            % Apply STDP for excitatory connections in the positive time case.
            % Uncomment this next line to disable STDP
           nuclei{receivingNuclei}{nucleiNumber}.S(:,excitatorySubset) = ...
               nuclei{receivingNuclei}{nucleiNumber}.S(:,excitatorySubset) ...
               + stdpEffectsPos(:,excitatorySubset);
            % Maximum weight size is an indication of the hypothetical minimum
            % number of incoming spikes required to cause a postsynaptic spike.
            % Clip weights to between upper bounds,
    		nuclei{receivingNuclei}{nucleiNumber}.S(nuclei{receivingNuclei}{nucleiNumber}.S ...
                > nuclei{receivingNuclei}{nucleiNumber}.params.weightUpperbound) ...
                = nuclei{receivingNuclei}{nucleiNumber}.params.weightUpperbound;
        end
        if ~isempty(nuclei{nucleiNumber}{receivingNuclei}.lastFire)
            % Apply STDP for excitatory connections in the negative time case.
            % % Uncomment this next line to disable STDP
           nuclei{nucleiNumber}{receivingNuclei}.S(excitatorySubset,:) = ...
               nuclei{nucleiNumber}{receivingNuclei}.S(excitatorySubset,:) ...
               + stdpEffectsNeg(excitatorySubset,:);
            % and lower bounds.
            nuclei{nucleiNumber}{receivingNuclei}.S(nuclei{nucleiNumber}{receivingNuclei}.S ...
                < nuclei{nucleiNumber}{receivingNuclei}.params.weightLowerbound) = ...
                nuclei{nucleiNumber}{receivingNuclei}.params.weightLowerbound;
        end

        if ~isempty(nuclei{nucleiNumber}{receivingNuclei}.lastFire)
            % Don't let excitatory synapses become inhibitory synapses by
            % STDP. This process does not need to be repeated for inhibitory
            % synapses because STDP should not be applied to inhibitory.
            % It does not need to applied to the nuclei{receivingNuclei}{nucleiNumber} because it is not receiving any negative STDP, so it could not have just become inhibitory.
            nuclei{nucleiNumber}{receivingNuclei}.S(nuclei{nucleiNumber}{receivingNuclei}.S(1:nuclei{nucleiNumber}{nucleiNumber}.params.Ne,1:nuclei{receivingNuclei}{receivingNuclei}.params.Ne) < 0) = 0;
        end

        % This is debug code for examining the case where STDP incorrectly brings a NaN in from lastFire.
		%       if nnz(isnan(stdpEffects)) > 0
		%           keyboard
		%       end
	end
end
clear receivingNuclei
%x(2) = toc;
%tic
% This is code remaining from a gating value implementation. It would need modification to be used, so it remains to make it easier to implement gating values in the future. This is part of the model conditions required to exhibit the izhi-szat working memory conditions. 
% Apply the STDP effects to the short-term STDP mechanism.
%          gatingValues(:,excitatorySubset) = gatingValues(:,excitatorySubset) ...
%              + stdpEffects(:,excitatorySubset)*100;
% ASSUMPTION: Gating values should not go below zero.
%          gatingValues(gatingValues < 0) = 0;



if ~isempty(fired)
  % Every spike causes N potentials. Each iteration of this for
  % loop adds N spike seperate potentials to conductingPotentials.
	for k = 1:length(fired)
		for receivingNuclei = 1:numNuclei
			% When a neuron fires, all of its dendrites immediately know, which is relevant for STDP.
			% Hence, we must update all of their lastFires.
			% We loop through each possible sending nuclei.
			% Only update lastFires where there actually exist connections to the nuclei of interest! (nucleiNumber)
			% For this statement, we use receivingNuclei as if it is the "sendingNuclei", so the self-documenting variable name is misleading for this statement.
			if ~isempty(nuclei{receivingNuclei}{nucleiNumber}.lastFire)
				% There exist connections, so update the relevant lastFire entries.
				nuclei{receivingNuclei}{nucleiNumber}.lastFire(nuclei{receivingNuclei}{nucleiNumber}.connectivityMatrix(:,fired(k)) == 1,fired(k)) = t;
            end
			
            if ~isempty(nuclei{nucleiNumber}{receivingNuclei}.lastFire)
                % Determine the relevant outgoing connections.
                effectiveConnections = find(nuclei{nucleiNumber}{receivingNuclei}.connectivityMatrix(fired(k),:));
                % Assemble a list of events for all nonzero outgoing
                % connections.
                newPotentials = [t*ones(length(effectiveConnections),1)  ... %the current time + how long it takes to send is
                  + nuclei{nucleiNumber}{receivingNuclei}.conductanceDelays(fired(k),effectiveConnections)', ... %when they will arrive
                ones(length(effectiveConnections),1)*fired(k), ... %where they came from (neuron number)
                transpose(effectiveConnections)]; % where they are going (neuron number)
            
            % This is a debug test...
            % if the new potentials should arrive at the current time,
            % something is wrong.
   %         if nnz(newPotentials(:,1) == t)>0
   %             keyboard
   %        end

                % Combine the new list with the old list.
                nuclei{nucleiNumber}{receivingNuclei}.conductingPotentials = [nuclei{nucleiNumber}{receivingNuclei}.conductingPotentials; newPotentials];
            end
		end
	end
end
%x(3) = toc;
%tic

% If no potentials have arrived, input will stay zero. Otherwise,
% it will grow.
input = zeros(nuclei{nucleiNumber}{nucleiNumber}.params.N,1);
% Find all received potentials for each sending nuclei
for sendingNuclei = 1:numNuclei
    if ~isempty(nuclei{sendingNuclei}{nucleiNumber}.lastFire)
        if ~isempty(nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials)
          % row indices of all potentials that have arrived at their
          % destination.
          arrivedPotentials = nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials(:,1) == t;
          if ~isempty(arrivedPotentials)
              % We create a potentialList to minimize the size of the index.
              % Since conductingPotentials can grow quite large, we want to
              % only reference it the minimum number of times. Particularly, I
              % have noticed that using syntax like
              % conductingPotentials(arrivedPotentials,:) is very
              % computationally expensive and has enhanced huge slow downs.
              potentialList = nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials(arrivedPotentials,:);
              % for each arrived potential
              for potIndex = 1:size(potentialList,1)
                  % The potential is being received by a neuron, so update the lastFire term.
                  nuclei{sendingNuclei}{nucleiNumber}.lastFire(potentialList(potIndex,2), ...
                           potentialList(potIndex,3)) = t;

                  % Calculate the input from the potential.
                  % In izhi&szat, I = S*(1+sd), where sd is a term that
                  % grows and shrinks by a form of STDP, and has constant
                  % decay.
                  % without gatingValues
                  input(potentialList(potIndex,3)) = ...
                      + input(potentialList(potIndex,3)) ...
                      + nuclei{sendingNuclei}{nucleiNumber}.S(potentialList(potIndex,2), potentialList(potIndex,3));

                  % with gatingValues
        %                       input(potentialList(potIndex,3)) = ...
        %                           + input(potentialList(potIndex,3)) ...
        %                           + S(potentialList(potIndex,2), potentialList(potIndex,3)) ...
        %                           .* (1 + gatingValues(potentialList(potIndex,2), potentialList(potIndex,3)));
              end
                % Clip out arrived action potentials.
                nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials = nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials(~arrivedPotentials, :);
          end
        end
    end
end
%x(4) = toc;
%tic

% inputList is a planned sequence of firings such that the per second
% firing rate of the 1000 neurons is 0.3, as implied by izhi&szat.
%synapticNoise = inputList(inputList(:,1) == t, 2);

inputNoise = zeros(nuclei{nucleiNumber}{nucleiNumber}.params.N,1);
inputNoise(willFire) = noiseAmplitude;

%synapticNoise = rand(N,1)*4;
% include thalamic input
input=input + inputNoise;

% Membrane Potential Update.
% step 0.5 ms for numerical stability
%  v=v+0.5*((0.04*v+5).*v+140-u+I);
nuclei{nucleiNumber}{nucleiNumber}.membranePotential = nuclei{nucleiNumber}{nucleiNumber}.membranePotential + nuclei{nucleiNumber}{nucleiNumber}.params.timeStep * ((0.04*nuclei{nucleiNumber}{nucleiNumber}.membranePotential+5).*nuclei{nucleiNumber}{nucleiNumber}.membranePotential + 140 - nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce + input);
nuclei{nucleiNumber}{nucleiNumber}.membranePotential = nuclei{nucleiNumber}{nucleiNumber}.membranePotential + nuclei{nucleiNumber}{nucleiNumber}.params.timeStep * ((0.04*nuclei{nucleiNumber}{nucleiNumber}.membranePotential+5).*nuclei{nucleiNumber}{nucleiNumber}.membranePotential + 140 - nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce + input);
% u=u+a.*(0.2*v-u);
nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce = nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce + nuclei{nucleiNumber}{nucleiNumber}.params.a.* (nuclei{nucleiNumber}{nucleiNumber}.params.b.*nuclei{nucleiNumber}{nucleiNumber}.membranePotential - nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce);

% This accomplishes an exponential decay which nears 0.0183 (close to
% zero) after 5000 iterations. 0.0183 = exp(5000 * log(.9992)). I
% chose the factor .9992 visually after plotting many points, so as
% to best fit the description in izhi&szat: "it decays back to 0 with
% a time constant 5 seconds." 
%          gatingValues = gatingValues*.9992;
% Gating Values have a maximum of 100% increase of excitation.
% "A maximum of 100% temporary increase relative to baseline"
%          gatingValues(gatingValues > 1) = 1;


%x(5) = toc;

%clear arrivedPotentials effectiveConnections excitatorySubset fired
%clear k newPotentials potIndex potentialList sendingNuclei
%clear stdpEffectsNeg stdpEffectsPos input