function [ newDataset ] = trimDataset( dataset , nrOfSamples )
%trimDataset Trim a dataset of the first and last nrOfSamples samples
    newDataset = dataset;
    
    % save timestamps
    newDataset.timeStamp = dataset.timeStamp(nrOfSamples:end-nrOfSamples);
    
    % save joint positions
    newDataset.q = dataset.q(nrOfSamples:end-nrOfSamples,:);
    
    % save joint velocities
    newDataset.dq = dataset.dq(nrOfSamples:end-nrOfSamples,:);
    
    % save joint accelerations 
    newDataset.ddq = dataset.ddq(nrOfSamples:end-nrOfSamples,:);
    
    % save sensor measures 
    %newDataset.ft = dataset.ft(nrOfSamples:end-nrOfSamples,:);
    
    % save sensor measures without termal drift
    %newDataset.ftNoDrift = dataset.ftNoDrift(nrOfSamples:end-nrOfSamples,:);
end

