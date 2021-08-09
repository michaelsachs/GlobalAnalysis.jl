using PlotlyJS
using Dash, DashHtmlComponents, DashCoreComponents, DashTable, DashDaq

# access at http://127.0.0.1:8050/

# note that Dash uses 0-indexing, hence some converting to/from Julia indexes is required

# issue with most recently PlotlyJS version: Failed component prop type: Invalid component prop `figure` key `config` supplied to Graph.
# using v0.14.0 instead


#numDat = 2

function normTrc(trc)
    maxi = abs.(maximum(trc, dims=1))
    mini = abs.(minimum(trc, dims=1))
    trc ./ (maxi > mini ? maxi : mini)
end

#cell type is NamedTuple
trc1(cell) = PlotlyJS.contour(; 
    x = Data[cell.row+1].x, 
    y = Data[cell.row+1].y, 
    #transpose for standalone PlotlyJS, but not for Dash 
    z = Data[cell.row+1].z,
    name = "data", 
    colorscale="Viridis", 
    line_width=0,
    ncontours = 40)
    
lay1 = Layout(;
    title= "Data",
    xaxis_title= "Wavelength (nm)",
    yaxis_title= "Time (ps)",
    yaxis_type= "log",
    grid_rows = 1,
    grid_columns = 1,
    height = 800
)
#PlotlyJS.plot(trc1(Dict("row" => 0)), lay1)


trc2(ystr,numDat,normBool) = vcat(
    [map(idx -> PlotlyJS.scatter(; 
        x = Data[n].x, 
        y = normBool ? normTrc(Data[n].z[:,idx]) : Data[n].z[:,idx],
        name = Data[n].y[idx], 
        colorscale="Viridis", 
        line_width=2), 
    findIdx(parseInputString(ystr),Data[n].y)) for n in numDat]...)

lay2 = Layout(;
    title= "Spectra",
    xaxis_title= "Wavelength (nm)",
    yaxis_title= "ΔA",
    xaxis_showline = true,
    xaxis_ticks = "outside",
    xaxis_mirror = true,
    xaxis_gridcolor = "#E8E8E8",
    yaxis_showline = true,
    yaxis_ticks = "outside",
    yaxis_mirror = true,
    yaxis_gridcolor = "#E8E8E8",
    height = 400
    )

trc3(xstr,numDat,normBool) = vcat(
    [map(idx -> PlotlyJS.scatter(; 
        x = Data[n].y, 
        y = normBool ? normTrc(Data[n].z[idx,:]) : Data[n].z[idx,:],
        name = Data[n].x[idx], 
        colorscale="Viridis", 
        line_width=2), 
    findIdx(parseInputString(xstr),Data[n].x)) for n in numDat]...)

lay3 = Layout(;
    title= "Kinetics",
    xaxis_title= "Time (ps)",
    yaxis_title= "ΔA",
    xaxis_type= "log",
    xaxis_showline = true,
    xaxis_ticks = "outside",
    xaxis_mirror = true,
    xaxis_gridcolor = "#E8E8E8",
    yaxis_showline = true,
    yaxis_ticks = "outside",
    yaxis_mirror = true,
    yaxis_gridcolor = "#E8E8E8",
    height = 400
    )
#Plot(trc3a([80,90]),lay3)

function parseInputString(str)
    parse.(Float64, split(str, " "))
end

# find closest indices of elements of val in vector vec
function findIdx(val, vec)
    idc = Int64[]
    for n in 1:length(val)
        idx = findmin(abs.(vec .- val[n]))[2]
        push!(idc, idx)
    end
    return idc
end


app = dash()

app.layout = html_div() do
    #html_h4("Explore data"),

    #html_button("Light", id="btn_theme", n_clicks=0),


    #dcc_dropdown(id = "dd_theme",
    #    options = [
    #        (label = "Light", value = "theme_light"),
    #        (label = "Dark", value = "theme_dark")
    #    ],
    #    value = "theme_light",
    #),

    html_div(id="table1", style=(backgroundColor="#dbe4f0 ", border="4px solid LightSteelBlue",hidden=false), #width=default_width,
        DashTable.dash_datatable(id="table", editable=true,
            columns = [(name="name", id="name"),(name="var", id="var")],
            data = [Dict("name" => Data[n].name, "var" => Data[n].var) for n in eachindex(Data)],
            row_selectable="multi",
            selected_rows=[0],
            active_cell = Dict("row" => 0, "column" => 0)
        )
    ),

    dcc_graph(id = "contour",
        #figure = Plot(trc1(1), lay1),
        style = (width = "49%", display = "inline-block")
    ),

    html_div(id="kinspc", style = (width = "49%", display = "inline-block"),
        children = [
            dcc_graph(
                id = "fig_spc",
                #figure = Plot(trc2, lay2),
                style = (width = "100%", display = "block")
            ),

            dcc_graph(
                id = "fig_kin",
                style = (width = "100%", display = "block")
            )
        ]
    ),

    dcc_input(id="input-spc", 
        type="text", 
        value="0.5 1 3 10 100 1000 6000",
        debounce=true),

    dcc_input(id="input-kin", 
        type="text", 
        value="580 700",
        debounce=true),

    daq_booleanswitch(id="norm_switch",
        on=false)
end


callback!(app,
    Output("contour", "figure"),
    Input("table", "active_cell")
) do active_cell
    return Plot(
        trc1(active_cell),
            lay1
        )
end


callback!(app,
    Output("fig_spc", "figure"),
    Input("input-spc", "value"),
    Input("table", "selected_rows"),
    Input("norm_switch", "on")
) do input, selected_rows, normBool
    return Plot(
        trc2(input, selected_rows .+ 1, normBool), 
            lay2
        )
end

callback!(app,
    Output("fig_kin", "figure"),
    Input("input-kin", "value"),
    Input("table", "selected_rows"),
    Input("norm_switch", "on")
) do input, selected_rows, normBool
    return Plot(
        trc3(input, selected_rows .+ 1, normBool), 
            lay3
        )
end


run_server(app, "0.0.0.0", debug=true)



