% parameter - (matlab.graphics.axis.Axes, double, double, double)
function func_plot_graph(axes_name, mat_cv, mat_rt, mat_intensity,bone)
    axes(axes_name);
    surf(mat_cv,mat_rt,mat_intensity);
    shading interp
    view(0,90)
    xlim([min(mat_cv),max(mat_cv)])
    ylim([min(mat_rt),max(mat_rt)])
    if nargin == 5
        colormap(axes_name,bone);
    end
end