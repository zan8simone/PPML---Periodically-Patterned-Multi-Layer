% plot_tiled_reflectance
% Make a tiled figure of 2D reflectance vs frequency plots
% for various angles of incidence

function plot_tiled_reflectance(R, angle, frequency)
    figure
    t = tiledlayout("horizontal");
    title(t, "Reflectance for incidence angles")
    xlabel(t, "Reflectance")
    ylabel(t, "Frequency (THz)")
    for angle_index = 1:length(angle)
        nexttile
        plot(R(:, angle_index), frequency)
        title([num2str(angle(angle_index)), char(176)])
    end

end