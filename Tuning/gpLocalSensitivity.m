function [gradMu, gradSigma2] = gpLocalSensitivity(gp, x_eval)
% Calculates the exact gradient of the GP predictive mean and variance at x_eval (original scale)
% Works with fitrgp using 'ardsquaredexponential' kernel

    % 1. Standardize input if necessary
    if isprop(gp, 'InputMu') && isprop(gp, 'InputSigma')
        x_std = (x_eval - gp.InputMu) ./ gp.InputSigma;
    else
        x_std = x_eval;
    end

    % 2. Extract training data and parameters (all standardized)
    X = table2array(gp.X);
    y = gp.Y;
    D = numel(x_std);
    N = size(X,1);

    kernelParams = gp.KernelInformation.KernelParameters;
    lengthscales = kernelParams(1:D).';
    sigmaF = kernelParams(end);
    sigmaN = gp.Sigma;

    % Target mean (center)
    if isprop(gp, 'Mu')
        mu_y = gp.Mu;
    else
        mu_y = mean(y);
    end

    % 3. Kernel matrix K (N x N)
    K = zeros(N,N);
    for i = 1:N
        for j = 1:N
            diff = (X(i,:) - X(j,:)) ./ lengthscales;
            K(i,j) = sigmaF^2 * exp(-0.5 * sum(diff.^2));
        end
    end
    K = K + sigmaN^2 * eye(N);

    % 4. k_star vector (1 x N)
    k_star = zeros(1,N);
    for n = 1:N
        diff = (x_std - X(n,:)) ./ lengthscales;
        k_star(n) = sigmaF^2 * exp(-0.5 * sum(diff.^2));
    end

    % 5. alpha (N x 1)
    alpha = K \ (y - mu_y);

    % 6. Predictive mean
    mu_pred = k_star * alpha + mu_y;

    % 7. Predictive variance
    v = K \ k_star';
    sigma2_pred = sigmaF^2 - k_star * v;

    % 8. Gradient of mean (in standardized space)
    gradMu = zeros(1,D);
    for d = 1:D
        deriv_k = zeros(1,N);
        for n = 1:N
            deriv_k(n) = k_star(n) * (X(n,d) - x_std(d)) / (lengthscales(d)^2);
        end
        gradMu(d) = deriv_k * alpha;
    end

    % 9. Gradient of variance (in standardized space)
    gradSigma2 = zeros(1,D);
    for d = 1:D
        % Derivative of k_star wrt x_d
        deriv_k = zeros(1,N);
        for n = 1:N
            deriv_k(n) = k_star(n) * (X(n,d) - x_std(d)) / (lengthscales(d)^2);
        end
        % gradSigma2(d) = -2 * deriv_k * v; (ver fórmula más abajo)
        gradSigma2(d) = -2 * (deriv_k * v);
    end

    % 10. Convert gradients to original input scale (chain rule)
    if isprop(gp, 'InputSigma')
        gradMu     = gradMu     ./ gp.InputSigma;
        gradSigma2 = gradSigma2 ./ gp.InputSigma;
    end
end