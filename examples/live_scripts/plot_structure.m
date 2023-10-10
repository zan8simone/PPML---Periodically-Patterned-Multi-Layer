function plot_structure(f, d, a)
    % plot_structure Plot structure of material for visual check
    % f: fraction of material B in material A
    % d: depth of layer (microns)
    % a: stripe spacing (microns)
    figure
    set(gca, 'YDir','reverse')
    end_depth = cumsum(d);
    start_depth = [0, end_depth(1:end -1)];
    z_offset = start_depth(end);
    start_depth = -(start_depth - z_offset);
    end_depth = -(end_depth - z_offset);
    colours = ["b", "r", "g", "c", "m", "y", "k"];
    internal_layer_num = 0;
    for layer = 1:length(d)
        if layer == 1
            layer_name = "Substrate";
        elseif layer == length(d)
            layer_name = "Superstrate";
        else
            layer_name = "Internal layer";
        end
        yregion(start_depth(layer), end_depth(layer), FaceColor=colours(layer), DisplayName=layer_name, FaceAlpha=0.7)

        if layer_name == "Internal layer"
            internal_layer_num = internal_layer_num + 1;
            fraction = f(internal_layer_num);
            half_x_width = 0.5 * fraction * a;
            x1 = - half_x_width;
            x2 = half_x_width;
            x = [x1, x1, x2, x2];
            y = [start_depth(layer), end_depth(layer), end_depth(layer), start_depth(layer)];
            patch(x, y, 'k', 'HandleVisibility', 'off', 'FaceAlpha', 0.5)
        end
    end
    substrate_display_depth = 10;
    ylim([start_depth(end) - substrate_display_depth, (end_depth(1) + substrate_display_depth)])
    xlim([-a/2, a/2])
    mu = char(181);
    ylabel("z coordinate (" + mu + "m)")
    xlabel("x coordinate (" + mu + "m)")
    legend('Direction','reverse')
    title("Sample structure")
end

