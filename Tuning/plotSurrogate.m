function stop = plotSurrogate(results, state)
    stop = false;
    if ~strcmp(state, 'iteration')
        return;
    end

    varsDesc = results.VariableDescriptions;
    varNames = {varsDesc.Name};
    D        = numel(varNames);
    Xobs     = results.XTrace;
    Yobs     = results.ObjectiveTrace;
    bestEst  = table2struct(results.XAtMinEstimatedObjective);

    % Fit the GP model
    gp = fitrgp(Xobs, Yobs, ...
        'Standardize', true, ...
        'KernelFunction', 'ardsquaredexponential');

    % ====== 1. UNIVARIATE SLICES (conditioned at best point) ======
    figure(5); clf;
    for iVar = 1:D
        lb = varsDesc(iVar).Range(1);
        ub = varsDesc(iVar).Range(2);
        xgrid = logspace(log10(lb), log10(ub), 120)';
        Tpred = repmat(struct2table(bestEst), numel(xgrid), 1);
        Tpred.(varNames{iVar}) = xgrid;
        [mu, sigma] = predict(gp, Tpred);

        subplot(D, 1, iVar);
        hold on; grid on;
        fill([xgrid; flipud(xgrid)], [mu+2*sigma; flipud(mu-2*sigma)], 0.9*[1 1 1], 'EdgeColor','none');
        plot(xgrid, mu, '-k', 'LineWidth', 1.5);
        scatter(Xobs{:,iVar}, Yobs, 30, 'r', 'filled');
        set(gca, 'XScale', 'log');
        xlabel(varNames{iVar}, 'Interpreter','none');
        ylabel('Cost');
    end
    sgtitle('1D GP slices, conditioned at best point', 'FontWeight', 'bold');

    % ==== 2. ONLY lambda_z vs lambda_o (conditioned on best lambda), SLICES UP/RIGHT ====
    figure(6); clf;
    idx_lambda_z = find(strcmp(varNames, 'lambda_z'));
    idx_lambda_o = find(strcmp(varNames, 'lambda_o'));
    idx_lambda   = find(strcmp(varNames, 'lambda'));

    % Grids for lambda_z and lambda_o
    lb_z = varsDesc(idx_lambda_z).Range(1); ub_z = varsDesc(idx_lambda_z).Range(2);
    lb_o = varsDesc(idx_lambda_o).Range(1); ub_o = varsDesc(idx_lambda_o).Range(2);
    xi   = logspace(log10(lb_z), log10(ub_z), 70);
    xj   = logspace(log10(lb_o), log10(ub_o), 70);
    [Xi, Xj] = meshgrid(xi, xj);

    % --- Layout: [left bottom width height] ---
    mainPos   = [0.13, 0.13, 0.60, 0.62];   % central (reduced width)
    topPos    = [0.13, 0.76, 0.60, 0.17];   % top (same width as main)
    rightPos  = [0.74, 0.13, 0.13, 0.62];   % right (move left, reduce width)
    cbPos     = [0.88, 0.13, 0.025, 0.62];  % colorbar (slim and separated)

    % --- Central bivariate contour (lambda_z vs lambda_o) ---
    Tcont = repmat(struct2table(bestEst), numel(Xi), 1);
    Tcont.(varNames{idx_lambda_z}) = Xi(:);
    Tcont.(varNames{idx_lambda_o}) = Xj(:);
    MuC = predict(gp, Tcont);
    MuC = reshape(MuC, size(Xi));

    axMain = axes('Position', mainPos);
    contourf(axMain, Xi, Xj, MuC, 20, 'LineColor','none');
    hold(axMain, 'on');
    scatter(axMain, Xobs{:,idx_lambda_z}, Xobs{:,idx_lambda_o}, 25, Yobs, 'filled', 'MarkerEdgeColor','k');
    % --------- RESALTAR PUNTO OPTIMO EN ROJO -----------
    % Obtener el punto óptimo (mínimo costo)
    [~, idx_best] = min(Yobs);
    lambda_z_opt = Xobs{idx_best, idx_lambda_z};
    lambda_o_opt = Xobs{idx_best, idx_lambda_o};
    plot(axMain, lambda_z_opt, lambda_o_opt, 'ro', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    % ---------------------------------------------------
    set(axMain, 'XScale','log', 'YScale','log');
    xlabel(axMain, 'lambda\_z');
    ylabel(axMain, 'lambda\_o');
    grid(axMain, 'on');
    box(axMain, 'on');

    % --- Top: lambda_z slice ---
    xgrid_z = logspace(log10(lb_z), log10(ub_z), 100)';
    Tz = repmat(struct2table(bestEst), numel(xgrid_z), 1);
    Tz.(varNames{idx_lambda_z}) = xgrid_z;
    [mu_z, sigma_z] = predict(gp, Tz);

    axTop = axes('Position', topPos);
    fill(axTop, [xgrid_z; flipud(xgrid_z)], ...
         [mu_z+2*sigma_z; flipud(mu_z-2*sigma_z)], 0.9*[1 1 1], 'EdgeColor', 'none');
    hold(axTop, 'on');
    plot(axTop, xgrid_z, mu_z, '-k', 'LineWidth', 1.5);
    scatter(axTop, Xobs{:,idx_lambda_z}, Yobs, 25, 'r', 'filled');
    set(axTop, 'XScale','log');
    axTop.XAxisLocation = 'bottom';
    axTop.YAxisLocation = 'left';
    axTop.XTickLabel = {};
    axTop.YTickLabel = {};
    axTop.Box = 'on';
    grid(axTop, 'on');
    ylabel(axTop, '');

    % --- Right: lambda_o slice ---
    xgrid_o = logspace(log10(lb_o), log10(ub_o), 100)';
    To = repmat(struct2table(bestEst), numel(xgrid_o), 1);
    To.(varNames{idx_lambda_o}) = xgrid_o;
    [mu_o, sigma_o] = predict(gp, To);

    axRight = axes('Position', rightPos);
    fill(axRight, [mu_o+2*sigma_o; flipud(mu_o-2*sigma_o)], ...
         [xgrid_o; flipud(xgrid_o)], 0.9*[1 1 1], 'EdgeColor', 'none');
    hold(axRight, 'on');
    plot(axRight, mu_o, xgrid_o, '-k', 'LineWidth', 1.5);
    scatter(axRight, Yobs, Xobs{:,idx_lambda_o}, 25, 'r', 'filled');
    set(axRight, 'YScale','log');
    axRight.XAxisLocation = 'top';
    axRight.YAxisLocation = 'right';
    axRight.XTickLabel = {};
    axRight.YTickLabel = {};
    axRight.Box = 'on';
    grid(axRight, 'on');
    xlabel(axRight, '');

    % --- Colorbar ---
    axes(axMain);
    c = colorbar('Position', cbPos);
    c.Label.String = 'Cost';

    sgtitle('2D GP slice, conditioned at best point', 'FontWeight', 'bold');
    drawnow;
end
