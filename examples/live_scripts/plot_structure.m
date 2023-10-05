function plot_structure(L, f, d, a)
    % plot_structure Plot structure of material for visual check
    % patch to fill between lines
    % or yregion
    figure
    end_depth = cumsum(d);
    start_depth = [0, end_depth(1:end -1)];
    colours = ["r", "b", "g", "c", "m", "y", "k"];
    for layer = 1:length(d)
        if layer == 1
            layer_name = "Substrate";
        elseif layer == length(d)
            layer_name = "Superstrate";
        else
            layer_name = "Internal layer";
        end

        yregion(start_depth(layer), end_depth(layer), FaceColor=colours(layer), DisplayName=layer_name)
    end
    ylim([start_depth(1), end_depth(end)])
    xlim([0, a])

    ylabel("Layer depth")
    legend()
end

