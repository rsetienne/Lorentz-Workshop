#julia --startup-file=no

using Unitful
using CarboKitten
using CairoMakie

facies = [
  ALCAP.Facies(
    viability_range=(4, 10),
    activation_range=(6, 10),
    production=BenthicProduction(
      maximum_growth_rate=500u"m/Myr",
      extinction_coefficient=0.8u"m^-1",
      saturation_intensity=60u"W/m^2"),
    transport_coefficient=50.0u"m/yr"),
  ALCAP.Facies(
    viability_range=(4, 10),
    activation_range=(6, 10),
    production=BenthicProduction(
      maximum_growth_rate=400u"m/Myr",
      extinction_coefficient=0.1u"m^-1",
      saturation_intensity=60u"W/m^2"),
    transport_coefficient=25.0u"m/yr"),
  ALCAP.Facies(
    viability_range=(4, 10),
    activation_range=(6, 10),
    production=BenthicProduction(
      maximum_growth_rate=100u"m/Myr",
      extinction_coefficient=0.005u"m^-1",
      saturation_intensity=60u"W/m^2"),
    transport_coefficient=12.5u"m/yr")
]

box = CarboKitten.Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m")
time = TimeProperties(Δt=0.0002u"Myr", steps=5000)
initial_topography = (x, y) -> -x / 300.0

locations = [
  (key = :col_proximal, label = "proximal ramp",   x = 10, y = 25),
  (key = :col_midplat,  label = "mid-platform",    x = 35, y = 25),
  (key = :col_margin,   label = "platform margin", x = 65, y = 25),
  (key = :col_slope,    label = "distal slope",    x = 90, y = 25),
]

output = Dict(
  :topography => OutputSpec(write_interval=100),
  :profile    => OutputSpec(slice=(:, 25)),
)

for loc in locations
   output[loc.key] = OutputSpec(slice=(loc.x, loc.y), write_interval=1)
end

sea_level = t -> 5.5u"m" * sin(2π * t / 0.18u"Myr")
subsidence_rate = 38.0u"m/Myr"
disintegration_rate = 80.0u"m/Myr"

input = ALCAP.Input(
  box = box,
  ca_interval = 1,
  time = time,
  output = output,
  initial_topography = initial_topography,
  sea_level = sea_level,
  subsidence_rate = subsidence_rate,
  insolation = 400.0u"W/m^2",
  facies = facies,
  depositional_resolution = 0.5u"m",
  sediment_buffer_size = 50,
  disintegration_rate = disintegration_rate,
  lithification_time = 100u"yr",
)

FILENAME = "data/output/carbo-platform.h5"
#using Downloads
#mkpath(dirname(FILENAME))
#Downloads.download(
#  "https://zenodo.org/records/22893028/files/carbo-platform.h5?download=1",
#  FILENAME)

mkpath("data/output")
run_model(Model{ALCAP}, input, FILENAME)

using CarboKitten.Export: read_slice
using CarboKitten.Visualization: sediment_profile

header, section = read_slice(FILENAME, :profile)
sediment_profile(header, section)

using CarboKitten.Export: read_column, extract_sac, extract_sc, age_depth_model

function column_bands(filename, loc)
    header, data = read_column(filename, loc.key)
    sac_df = extract_sac(header, data, loc.key)
    adm = age_depth_model(sac_df)[:, "adm_$(loc.key)"] |> in_units_of(u"m")

    sc_df = extract_sc(header, data, loc.key)
    facies_cols = ["sc_$(loc.key)_f$(f)" for f in 1:header.n_facies]
    sc = Matrix(sc_df[:, facies_cols])
    dominant = getindex.(argmax(sc, dims=2)[:, 1], 2)

    return adm, dominant
end

colors = Makie.wong_colors()
fig_logs = Figure(size=(700, 500))
axs = Axis[]
for (i, loc) in enumerate(locations)
    adm, dominant = column_bands(FILENAME, loc)
    ax = Axis(fig_logs[1, i], title=loc.label, xticksvisible=false, xticklabelsvisible=false)
    hspan!(ax, adm[1:end-1], adm[2:end]; color=colors[dominant[1:end-1]])
    push!(axs, ax)
end
axs[1].ylabel = "stratigraphic height [m]"
linkyaxes!(axs...)

Legend(fig_logs[2, 1:4],
    [PolyElement(color=colors[f]) for f in 1:3],
    ["euphotic (shallow, fast)", "oligophotic (mid depth)", "aphotic (deep, slow)"],
    orientation=:horizontal, tellwidth=false)

fig_logs

fig_adm = Figure(size=(700, 450))
ax = Axis(fig_adm[1, 1], xlabel="time [Myr]", ylabel="stratigraphic height [m]",
          title="Age-depth models across the platform")

for loc in locations
    header, data = read_column(FILENAME, loc.key)
    t = header.axes.t |> in_units_of(u"Myr")
    sac_df = extract_sac(header, data, loc.key)
    adm = age_depth_model(sac_df)[:, "adm_$(loc.key)"] |> in_units_of(u"m")
    lines!(ax, t, adm, label=loc.label)
end

axislegend(ax, position=:rb)
fig_adm

using CarboKitten.Export: CSV, data_export

cols = Dict(loc.key => read_column(FILENAME, loc.key)[2] for loc in locations)
header, _ = read_column(FILENAME, first(locations).key)

spec = CSV(
    :sediment_accumulation_curve => "data/output/columns_sac.csv",
    :age_depth_model             => "data/output/columns_adm.csv",
    :stratigraphic_column        => "data/output/columns_sc.csv",
    :water_depth                 => "data/output/columns_wd.csv",
    :metadata                    => "data/output/columns_metadata.toml")

data_export(spec, header, cols)
