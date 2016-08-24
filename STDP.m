function [posDS negDS] = STDP(lastFire1, fired, t, Ne2, lastFire2)
% function [posDS negDS] = STDP(lastFire1, fired, t, Ne2, lastFire2)
% lastFire1 : is a matrix that contains the last time each neuron
%    experienced each other neuron fire.
% fired : is a logical index for all of the neurons that just fired.
% t : is the current time.
% Ne2: number of excitatory neurons in the network that the columns of
% lastFire1 refer to.
% lastFire2 is the lastFire matrix where Ne2 refers to the presynaptic neurons and Ne1 refers to the postsynaptic neurons.
% posDS and negDS : the return variables.
% posDS should be the size of lastFire2, and negDS should be the size of lastFire1. 

% If a neuron has never fired, it has lastFire value of NaN. However, since
% NaN elements never return true on numerical comparisons, the address
% times > 0 and times <= 0 is false for all NaN elements, hence no NaNs
% will end up in dS. NaN values get pruned out elsewhere.

% calculate the time difference between last spike arrival and current
% time.
% "Excitatory to inhibitory and all inhibitory connections are non-plastic"
% this is why we only calculate for neuron indexes 1:Ne.
%
% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012

if ~isempty(lastFire2)
	posTimes = t - lastFire2(1:Ne2,fired); % only for A_+.
else
	posTimes = [];
end

if ~isempty(lastFire1)
	negTimes = lastFire1(fired,1:Ne2) - t; % only for A_-.
else
	negTimes = [];
end

%times =  t - lastFire(:,fired);

%%%%% Currently: STDP only works for scenarios where postsynaptic neurons
%%%%% spike after the presynaptic neuron, and not before. This needs to be
%%%%% changed. The change I propose is to look at fired lastfires to
%%%%% determine if the last fire of the currently fired are
%%%%% anti-correlated. If it was, it gets negatively effected.

% initialize dS.
posDS = zeros(size(lastFire2));
negDS = zeros(size(lastFire1));

% Parameters taken from Izhikevich (2006)
tauMinus = 20; % 20 ms.
tauPlus = 20; % 20 ms.
Aminus = -0.12;
Aplus = 0.1;

% Approximation function taken from Izhikevich et. al (2004)
% when times == 0, it counts as being BEFORE, since minimum delay is 1 ms.
for fireNum = 1:length(fired)
    if nnz(posTimes>0) > 0
        posDS(posTimes(:,fireNum)>0,fired(fireNum)) = posDS(posTimes(:,fireNum)>0,fired(fireNum)) +...
            Aplus * exp(-posTimes(posTimes(:,fireNum)>0,fireNum) / tauPlus);
    end
    if nnz(negTimes<=0) > 0
        negDS(fired(fireNum),negTimes(fireNum,:)<=0) = negDS(fired(fireNum),negTimes(fireNum,:)<=0) +...
        Aminus * exp(negTimes(fireNum, negTimes(fireNum,:)<=0) / tauMinus);
    end
end

return
end