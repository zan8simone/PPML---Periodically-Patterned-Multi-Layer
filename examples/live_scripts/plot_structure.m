function plot_structure(L, f, d, a)
    % plot_structure Plot structure of material for visual check
    % L: number of layers
    % f: fraction of material B in material A
    % d: depth of layer (microns)
    % a: stripe spacing (microns)
    figure
    end_depth = cumsum(d);
    start_depth = [0, end_depth(1:end -1)];
    colours = ["b", "r", "g", "c", "m", "y", "k"];
    x_mid = 0.5 * a;
    internal_layer_num = 0;
    for layer = 1:length(d)
        if layer == 1
            layer_name = "Substrate";
        elseif layer == length(d)
            layer_name = "Superstrate";
        else
            layer_name = "Internal layer";
        end
        yregion(start_depth(layer), end_depth(layer), FaceColor=colours(layer), DisplayName=layer_name, FaceAlpha=0.8)

        if layer_name == "Internal layer"
            internal_layer_num = internal_layer_num + 1;
            fraction = f(internal_layer_num);
            x1 = x_mid - 0.5 * fraction * a;
            x2 = x_mid + 0.5 * fraction * a;
            x = [x1, x1, x2, x2];
            y = [start_depth(layer), end_depth(layer), end_depth(layer), start_depth(layer)];
            patch(x, y, 'k', 'HandleVisibility', 'off', 'FaceAlpha', 0.5)
        end
    end
    ylim([start_depth(1), end_depth(end)])
    xlim([0, a])
    mu = char(181);
    ylabel("Layer depth (" + mu + "m)")
    xlabel("Stripe spacing (" + mu + "m)")
    legend('Direction','reverse')
end

