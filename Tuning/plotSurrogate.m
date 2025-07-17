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

    % Fit the GP model
    gp = fitrgp(Xobs, Yobs, ...
        'Standardize', true, ...
        'KernelFunction', 'ardsquaredexponential');

    % ===== 1. UNIVARIATE SLICES =====
    figure(5); clf;
    for iVar = 1:D
        lb = varsDesc(iVar).Range(1);
        ub = varsDesc(iVar).Range(2);
        xgrid = logspace(log10(lb), log10(ub), 120)';
        % Para 2 variables, solo necesitas la otra fija al valor medio
        otherIdx = 3 - iVar; % 1->2, 2->1
        otherVal = mean(Xobs{:,otherIdx});
        Tpred = array2table([zeros(numel(xgrid),2)], 'VariableNames', varNames);
        Tpred{:,iVar} = xgrid;
        Tpred{:,otherIdx} = otherVal;
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
    sgtitle('Univariate GP slices', 'FontWeight', 'bold');

    % ===== 2. BIVARIATE GP (contour) =====
    figure(6); clf;
    lb1 = varsDesc(1).Range(1); ub1 = varsDesc(1).Range(2);
    lb2 = varsDesc(2).Range(1); ub2 = varsDesc(2).Range(2);
    xi = logspace(log10(lb1), log10(ub1), 70);
    xj = logspace(log10(lb2), log10(ub2), 70);
    [Xi, Xj] = meshgrid(xi, xj);

    Tcont = array2table([Xi(:), Xj(:)], 'VariableNames', varNames);
    MuC = predict(gp, Tcont);
    MuC = reshape(MuC, size(Xi));

    contourf(Xi, Xj, MuC, 20, 'LineColor','none');
    hold on;
    scatter(Xobs{:,1}, Xobs{:,2}, 30, Yobs, 'filled', 'MarkerEdgeColor','k');
    % Resalta el punto óptimo
    [~, idx_best] = min(Yobs);
    plot(Xobs{idx_best,1}, Xobs{idx_best,2}, 'ro', 'MarkerSize', 12, ...
        'MarkerFaceColor', 'r', 'LineWidth', 2);

    set(gca, 'XScale','log', 'YScale','log');
    xlabel(varNames{1}, 'Interpreter','none');
    ylabel(varNames{2}, 'Interpreter','none');
    colorbar;
    grid on; box on;
    sgtitle('2D GP surrogate', 'FontWeight', 'bold');
    drawnow;
end