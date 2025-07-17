function plotGP(gp, results)
    varsDesc = results.VariableDescriptions;
    varNames = {varsDesc.Name};
    Xobs     = results.XTrace;
    Yobs     = results.ObjectiveTrace;

    % Asume solo dos variables
    idx1 = 1;
    idx2 = 2;
    lb1 = varsDesc(idx1).Range(1); ub1 = varsDesc(idx1).Range(2);
    lb2 = varsDesc(idx2).Range(1); ub2 = varsDesc(idx2).Range(2);

    xi = logspace(log10(lb1), log10(ub1), 70);
    xj = logspace(log10(lb2), log10(ub2), 70);
    [Xi, Xj] = meshgrid(xi, xj);

    % Armar tabla para predicción
    Tpred = array2table([Xi(:), Xj(:)], 'VariableNames', varNames);
    MuC = predict(gp, Tpred);
    MuC = reshape(MuC, size(Xi));

    figure('Name', '2D GP Surrogate', 'NumberTitle', 'off', 'Color', 'w');
    contourf(Xi, Xj, MuC, 20, 'LineColor','none');
    hold on;
    scatter(Xobs{:,1}, Xobs{:,2}, 30, Yobs, 'filled', 'MarkerEdgeColor','k');
    % Punto óptimo en rojo
    [~, idx_best] = min(Yobs);
    plot(Xobs{idx_best,1}, Xobs{idx_best,2}, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    set(gca, 'XScale','log', 'YScale','log');
    xlabel(varNames{1}, 'Interpreter','none');
    ylabel(varNames{2}, 'Interpreter','none');
    colorbar;
    title('2D GP surrogate');
    grid on; box on;
end
