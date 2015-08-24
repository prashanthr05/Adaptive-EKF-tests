function [ dataset ] = smoothAndEstimateVelAcc( dataset )
%smoothAndEstimateVelAcc smooth joint position and estimate vel and acc
%   

for channel = 1:size(dataset.q,2);
    % compute timestamp from dataset 
    y = dataset.q(:,channel);
    ds = mean(diff(dataset.timeStamp));
    nrOfSamples = length(dataset.timeStamp);
    N = 2;
    F = 15; % 1 second window
    [b,g] = sgolay(N,F);
    HalfWin  = ((F+1)/2) -1;

    for n = (F+1)/2:nrOfSamples-(F+1)/2,
        % Zeroth derivative (smoothing only)
        dataset.q(n,channel) = dot(g(:,1),y(n - HalfWin:n + HalfWin));

        % 1st differential
        dataset.dq(n,channel) = dot(g(:,2),y(n - HalfWin:n + HalfWin));

        % 2nd differential
        dataset.ddq(n,channel) = 2*dot(g(:,3)',y(n - HalfWin:n + HalfWin))';
    end

    dataset.dq(:,channel) = dataset.dq(:,channel)/ds;         % Turn differential into derivative
    dataset.ddq(:,channel) = dataset.ddq(:,channel)/(ds*ds);    % and into 2nd derivative
    
    % trim the dataset

end

end

