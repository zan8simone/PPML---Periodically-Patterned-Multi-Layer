function plot_structure(f, d, a)
    % plot_structure Plot structure of material for visual check
    % f: fraction of material B in material A
    % d: depth of layer (microns)
    % a: stripe spacing (microns)
    % All these arrays are ordered from superstrate to substrate,
    figure

    % Flip order of layers so substrate is at the bottom of the plot
    d = flip(d);
    f = flip(f);

    end_depth = cumsum(d);
    start_depth = [0, end_depth(1:end -1)];

    % Reverse z coordinate axis direction
    set(gca, 'YDir','reverse')

    % Offset z values so bottom of superstrate is z=0
    z_offset = start_depth(end);
    start_depth = -(start_depth - z_offset);
    end_depth = -(end_depth - z_offset);

    colours = ["b", "r", "g", "c", "m", "y", "k"];
    internal_layer_num = 0;
    substrate_display_depth = 10;
    ymin = start_depth(end) - substrate_display_depth;
    ymax = (end_depth(1) + substrate_display_depth);
    ylim([ymin,ymax])
    xlim([-a/2, a/2])

    for layer = 1:length(d)
        if layer == 1
            layer_name = "Substrate";
            text(0.45 * a/2, ymax - 0.5 * substrate_display_depth, layer_name, 'FontSize',12);
        elseif layer == length(d)
            layer_name = "Superstrate";
            text(0.45 * a/2, ymin + 0.5 * substrate_display_depth, layer_name, 'FontSize',12);
        else
            internal_layer_num = internal_layer_num + 1;
            fraction = f(internal_layer_num);
            half_x_width = 0.5 * fraction * a;
            x1 = - half_x_width;
            x2 = half_x_width;
            x = [x1, x1, x2, x2];
            y = [start_depth(layer), end_depth(layer), end_depth(layer), start_depth(layer)];
            patch(x, y, 'k', 'HandleVisibility', 'off', 'FaceAlpha', 0.5)
        end
        yregion(start_depth(layer), end_depth(layer), FaceColor=colours(layer), FaceAlpha=0.7)
    end

    mu = char(181);
    ylabel("z coordinate (" + mu + "m)")
    xlabel("x coordinate (" + mu + "m)")
    title("Sample structure")
end

