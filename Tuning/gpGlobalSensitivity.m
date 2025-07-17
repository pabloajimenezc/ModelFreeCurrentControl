function lengthscales = gpGlobalSensitivity(gp)
% Extracts the lengthscales (ARD) from a fitrgp model
% Returns a vector (one per variable)

    kernelParams = gp.KernelInformation.KernelParameters;
    % In ARD Squared Exponential, the first D params are the lengthscales
    D = numel(gp.PredictorNames);
    lengthscales = kernelParams(1:D)';
end