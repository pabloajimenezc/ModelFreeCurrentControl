function plotGP(gp, results)
    varsDesc = results.VariableDescriptions;
    varNames = {varsDesc.Name};
    D        = numel(varNames);
    Xobs     = results.XTrace;
    Yobs     = results.ObjectiveTrace;

    % Indices
    idx_lambda_z = find(strcmp(varNames, 'lambda_z'));
    idx_lambda_o = find(strcmp(varNames, 'lambda_o'));
    idx_lambda   = find(strcmp(varNames, 'lambda'));

    % Lambda range
    lb_lambda = varsDesc(idx_lambda).Range(1);
    ub_lambda = varsDesc(idx_lambda).Range(2);
    sliderVals = logspace(log10(lb_lambda), log10(ub_lambda), 50);

    % Ajusta tamaño y espacios
    f = figure('Name', '2D GP Slice with Lambda Slider', 'NumberTitle', 'off');
    set(f, 'Color', 'w', 'Position', [100, 100, 580, 540]);
    mainPos   = [0.13, 0.25, 0.70, 0.70];   % Más alto y ancho
    cbPos     = [0.87, 0.25, 0.035, 0.70];  % Más pegada y alta
    sliderPos = [0.13 0.12 0.70 0.045];     % Debajo del gráfico principal
    txtPos    = [0.13 0.06 0.70 0.045];     % Texto debajo del slider

    % Main axes
    axMain = axes('Position', mainPos, 'Parent', f);

    % Slider y texto
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', numel(sliderVals), 'Value', 1, ...
        'SliderStep', [1/(numel(sliderVals)-1) 5/(numel(sliderVals)-1)], ...
        'Units', 'normalized', ...
        'Position', sliderPos);
    txt = uicontrol('Style','text', 'Units','normalized', ...
        'Position',txtPos, 'String','', ...
        'FontSize', 12, 'HorizontalAlignment', 'center', 'BackgroundColor', 'w');

    function updatePlot(idx)
        % Valor óptimo (primer fila de results ordenado)
        Xopt = results.XTrace(1, :);
        lambda_opt = Xopt.(varNames{idx_lambda});
        lambda_z_opt = Xopt.(varNames{idx_lambda_z});
        lambda_o_opt = Xopt.(varNames{idx_lambda_o});
    
        lambda_val = sliderVals(idx);
        set(txt, 'String', sprintf('lambda = %.2g', lambda_val));
    
        lb_z = varsDesc(idx_lambda_z).Range(1); ub_z = varsDesc(idx_lambda_z).Range(2);
        lb_o = varsDesc(idx_lambda_o).Range(1); ub_o = varsDesc(idx_lambda_o).Range(2);
        xi   = logspace(log10(lb_z), log10(ub_z), 70);
        xj   = logspace(log10(lb_o), log10(ub_o), 70);
        [Xi, Xj] = meshgrid(xi, xj);
    
        % Generar predicciones
        Tzero = array2table(zeros(1,D), 'VariableNames', varNames);
        Tpred = repmat(Tzero, numel(Xi), 1);
        Tpred.(varNames{idx_lambda_z}) = Xi(:);
        Tpred.(varNames{idx_lambda_o}) = Xj(:);
        Tpred.(varNames{idx_lambda})   = repmat(lambda_val, numel(Xi), 1);
    
        MuC = predict(gp, Tpred);
        MuC = reshape(MuC, size(Xi));
    
        % Limpiar solo el axis principal
        cla(axMain);
        contourf(axMain, Xi, Xj, MuC, 20, 'LineColor','none');
        hold(axMain, 'on');
        scatter(axMain, Xobs{:,idx_lambda_z}, Xobs{:,idx_lambda_o}, 25, Yobs, 'filled', 'MarkerEdgeColor','k');
        % Resalta el punto óptimo en rojo (más grande)
        plot(axMain, lambda_z_opt, lambda_o_opt, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    
        set(axMain, 'XScale','log', 'YScale','log');
        xlabel(axMain, 'lambda\_z');
        ylabel(axMain, 'lambda\_o');
        grid(axMain, 'on');
        box(axMain, 'on');
        c = colorbar(axMain, 'Position', cbPos);
        c.Label.String = 'Cost';
    
        % Título muestra valor óptimo
        title(axMain, sprintf('2D GP slice, optimal lambda = %.2g', lambda_opt), 'FontWeight', 'bold', 'FontSize', 14);
    end


    updatePlot(1);
    slider.Callback = @(src,~) updatePlot(round(src.Value));
end