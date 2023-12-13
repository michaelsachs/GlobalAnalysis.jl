using Plots

#----------------------------------------------------------------------
# Default plot layout
#----------------------------------------------------------------------

# default fonts
fo1 = font(9,"Helvetica")
fo2 = font(11,"Helvetica")

# set default plot parameters
default(linewidth=2, grid=false, frame=:box, fg_legend=:transparent,
    tickfont=fo1, guidefont=fo2, titlefont=fo2, legendfont=fo1,
    size=(500,350), legend=false, tick_direction=:out, dpi=300, 
    plot_titlefont=fo2)