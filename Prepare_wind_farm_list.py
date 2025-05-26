import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, shape
from shapely.ops import unary_union
from pyproj import Transformer
import shapefile  
import shapely
from pathlib import Path

# Import UK wind farms from your extraction module
from UKWindFarmsExtraction import UK_filtered_wind_farms

# Parameters
dist_from_shore_limit = 10  # km
min_separation_m = 1852*5 

# Norwegian wind farms (concept areas)
Nor_wind_farms = [
    {'name': 'Vestavind B', 'lat': 61.0000, 'lon': 3.5417, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Vestavind C', 'lat': 60.4000, 'lon': 3.7167, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Vestavind D', 'lat': 60.3583, 'lon': 4.4208, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Vestavind E', 'lat': 59.1524, 'lon': 3.8424, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Vestavind F (Utsira Nord)', 'lat': 59.3142, 'lon':  4.4856, 'installed_capacity': 1500, 'development_status': 'Utlyst', 'country': 'Norway'},
    {'name': 'Sørvest A', 'lat': 57.9429, 'lon': 3.5285, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Sørvest B', 'lat': 57.3720, 'lon': 3.4208, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Sørvest C', 'lat': 57.0301, 'lon': 3.9190, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Sørvest D', 'lat': 56.4410, 'lon': 3.7851, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Sørvest E', 'lat': 57.4722, 'lon': 4.7085, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'},
    {'name': 'Sørlige Nordsjø II', 'lat': 56.7375, 'lon': 5.0000, 'installed_capacity': 1500, 'development_status': 'Tildelt', 'country': 'Norway'},
    {'name': 'Sønnavind A', 'lat': 57.5298, 'lon': 7.5599, 'installed_capacity': 1000.0, 'development_status': 'Skal Utredes', 'country': 'Norway'}
]

# Farms found from: https://wind-analytics.esgian.com/
# Locations from: https://www.gem.wiki/Main_Page
duch_wind_farms = [
    {'name': 'IJmuiden Ver Alpha', 'lat': 52.7464, 'lon': 3.5694, 'installed_capacity': 2000.0, 'development_status': 'Pre-construction', 'country': 'Netherlands'},
    {'name': 'IJmuiden Ver Beta', 'lat': 52.9046, 'lon': 3.5827, 'installed_capacity': 2000.0, 'development_status': 'Pre-construction', 'country': 'Netherlands'},
    {'name': 'Hollandse Kust West VI', 'lat': 52.6932, 'lon': 3.8686, 'installed_capacity': 1050.0, 'development_status': 'Pre-construction', 'country': 'Netherlands'},
    {'name': 'Hollandse Kust Noord V', 'lat': 52.6424, 'lon': 4.3471, 'installed_capacity': 759.0, 'development_status': 'Operating', 'country': 'Netherlands'},
    {'name': 'Gemini', 'lat': 54.0470,  'lon': 5.9811, 'installed_capacity': 600, 'development_status': 'Operating', 'country': 'Netherlands'} #https://www.gem.wiki/Gemini_Offshore_Wind_Park
]

german_wind_farms = [
    {'name': 'Amrumbank West', 'lat': 54.5200, 'lon': 7.7080, 'installed_capacity': 302, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Amrumbank_West_wind_farm
    {'name': 'Kaskasi', 'lat': 54.4900,  'lon': 7.7000, 'installed_capacity': 342, 'development_status': 'Operating', 'country': 'Germany'},# https://www.gem.wiki/Kaskasi_wind_farm
    {'name': 'Nordsee Ost', 'lat': 54.4440,  'lon': 7.6820, 'installed_capacity': 332, 'development_status': 'Operating', 'country': 'Germany'},#https://www.gem.wiki/Nordsee_Ost_wind_farm
    {'name': 'Meerwind Süd Ost', 'lat': 54.4020,  'lon': 7.7070, 'installed_capacity': 288, 'development_status': 'Operating', 'country': 'Germany'},#https://www.gem.wiki/Meerwind_S%C3%BCd_Ost_wind_farm
    {'name': 'Global Tech', 'lat':54.5000, 'lon': 6.3580 , 'installed_capacity': 400, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Global_Tech_wind_farm
    {'name': 'Albatros', 'lat': 54.4868,  'lon': 6.2832, 'installed_capacity': 112, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Albatros_wind_farm
    {'name': 'Deutsche Bucht', 'lat': 54.3036,  'lon': 5.7894, 'installed_capacity': 260, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Deutsche_Bucht_wind_farm
    {'name': 'Veja Mate', 'lat': 54.3210, 'lon': 5.8600, 'installed_capacity': 402, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Veja_Mate_wind_farm
    {'name': 'BARD', 'lat': 54.3550, 'lon': 5.9800, 'installed_capacity': 400, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/BARD_Offshore_wind_farm
    {'name': 'Gode wind 1', 'lat': 54.0160, 'lon': 6.9830, 'installed_capacity': 330, 'development_status': 'Operating', 'country': 'Germany'},#https://www.gem.wiki/Gode_wind_farm
    {'name': 'Gode wind 2', 'lat': 54.0750, 'lon': 7.0070, 'installed_capacity': 252, 'development_status': 'Operating', 'country': 'Germany'},
    {'name': 'Gode wind 3', 'lat': 54.0370,  'lon': 7.1100, 'installed_capacity': 265, 'development_status': 'Operating', 'country': 'Germany'},
    {'name': 'Nordsee One', 'lat': 54.4440, 'lon': 7.6820, 'installed_capacity': 332, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Nordsee_One_wind_farm
    {'name': 'Merkur', 'lat': 54.0564, 'lon': 6.5536, 'installed_capacity': 396, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Merkur_wind_farm
    {'name': 'Trianel Borkum', 'lat': 54.0417, 'lon': 6.4667, 'installed_capacity': 402, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Trianel_Borkum_wind_farm
    {'name': 'Sandbank', 'lat': 55.1900, 'lon': 6.8600, 'installed_capacity': 288, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Sandbank_wind_farm
    {'name': 'DanTysk', 'lat': 55.1400,  'lon': 7.2000, 'installed_capacity': 288, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Dantysk_wind_farm
    {'name': 'Butendiek', 'lat': 55.0190,  'lon': 7.7740, 'installed_capacity': 240, 'development_status': 'Operating', 'country': 'Germany'}, #https://www.gem.wiki/Butendiek_wind_farm
    {'name': 'Borkum Riffgrund 1', 'lat': 53.9667, 'lon': 6.5623, 'installed_capacity': 420, 'development_status': 'Construction', 'country': 'Germany'}, #https://www.gem.wiki/Borkum_Riffgrund_3_wind_farm
    {'name': 'Borkum Riffgrund 2', 'lat': 53.9667, 'lon': 6.4956, 'installed_capacity': 240, 'development_status': 'Construction', 'country': 'Germany'},
    {'name': 'Borkum Riffgrund 3', 'lat': 54.0470, 'lon': 6.2340, 'installed_capacity': 240, 'development_status': 'Construction', 'country': 'Germany'}
    
]

danish_wind_farms = [
    {'name': 'Horns Rev 1', 'lat': 55.4691,  'lon': 7.8847, 'installed_capacity': 160, 'development_status': 'Operating', 'country': 'Denmark'}, #https://www.gem.wiki/Horns_Rev_Offshore_wind_farm
    {'name': 'Horns Rev 2', 'lat': 55.6084,   'lon': 7.6114, 'installed_capacity': 209, 'development_status': 'Operating', 'country': 'Denmark'},
    {'name': 'Horns Rev 3', 'lat': 55.6825,  'lon': 7.7729, 'installed_capacity': 407, 'development_status': 'Operating', 'country': 'Denmark'},
    {'name': 'Vesterhav Nord', 'lat': 56.6999, 'lon': 8.0446, 'installed_capacity': 180, 'development_status': 'Operating', 'country': 'Denmark'}, #https://www.gem.wiki/Vesterhav_Offshore_wind_farm
    {'name': 'Vesterhav Syd', 'lat': 56.0329,  'lon': 8.0267, 'installed_capacity': 170, 'development_status': 'Operating', 'country': 'Denmark'},
    {'name': 'Thor', 'lat': 56.3695, 'lon': 8.0143, 'installed_capacity': 1000, 'development_status': 'Pre-construction', 'country': 'Denmark'} #https://www.gem.wiki/Thor_Offshore_wind_farm
]

# Combine UK, Norwegian, and German wind farms
wind_farms = UK_filtered_wind_farms + Nor_wind_farms + german_wind_farms + danish_wind_farms + duch_wind_farms

# Convert the list of wind farms to a DataFrame
df = pd.DataFrame(wind_farms)
print(f"Number of wind farms before filtering: {len(df)}")

# Ensure latitude and longitude are numeric
df['lat'] = pd.to_numeric(df['lat'], errors='coerce')
df['lon'] = pd.to_numeric(df['lon'], errors='coerce')
df.dropna(subset=['lat', 'lon'], inplace=True)

# Create a geometry column from longitude and latitude
df['geometry'] = df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)

# Convert the DataFrame to a GeoDataFrame (CRS: WGS84)
gdf = gpd.GeoDataFrame(df, geometry='geometry', crs='EPSG:4326')

# Read the coastline shapefile
sf = shapefile.Reader("/Users/markusagersborg/Library/CloudStorage/OneDrive-NTNU/Master/Optimeringsmodell/ne_10m_coastline/ne_10m_coastline.shp")
coastline_shapes = sf.shapes()
coastline_geometries = [shape(s.__geo_interface__) for s in coastline_shapes]
coastline_union = unary_union(coastline_geometries)

# Set up coordinate transformations from WGS84 to a projected CRS (EPSG:32633)
transformer_to_proj = Transformer.from_crs("EPSG:4326", "EPSG:32633", always_xy=True)
transformer_to_wgs84 = Transformer.from_crs("EPSG:32633", "EPSG:4326", always_xy=True)

# Transform wind farm points to the projected CRS
def transform_point(point):
    x, y = transformer_to_proj.transform(point.x, point.y)
    return Point(x, y)

gdf['geometry_proj'] = gdf['geometry'].apply(transform_point)

# Transform coastline geometry to the projected CRS
coastline_union_proj = shapely.ops.transform(transformer_to_proj.transform, coastline_union)

# Calculate distance from each wind farm to the coastline (in meters)
def calculate_distance(point):
    return point.distance(coastline_union_proj)

gdf['distance_to_shore_m'] = gdf['geometry_proj'].apply(calculate_distance)
gdf['distance_to_shore_km'] = gdf['distance_to_shore_m'] / 1000

# Filter out wind farms that are closer than the specified distance from shore
filtered_gdf = gdf[gdf['distance_to_shore_km'] >= dist_from_shore_limit]
filtered_gdf = filtered_gdf[filtered_gdf['lon'] < 10]  # Filter out points east of 10 degrees longitude
print(f"Number of wind farms after shore filtering: {len(filtered_gdf)}")

pruned_indices = []
for idx, row in filtered_gdf.iterrows():
    current_point = row['geometry_proj']
    # Check against points already in pruned_indices
    if any(current_point.distance(filtered_gdf.loc[other_idx, 'geometry_proj']) < min_separation_m for other_idx in pruned_indices):
        continue  # Skip this point because it's too close to an already kept point
    pruned_indices.append(idx)

pruned_gdf = filtered_gdf.loc[pruned_indices]
print(f"Number of wind farms after proximity filtering (1 nm): {len(pruned_gdf)}")

# Convert back to a list of dictionaries (if needed)
wind_farm_list = pruned_gdf.drop(columns=['geometry', 'geometry_proj', 'distance_to_shore_m', 'distance_to_shore_km']).to_dict(orient='records')

# Calculate the center of the map (average latitude and longitude)
latitudes = [wf['lat'] for wf in wind_farm_list]
longitudes = [wf['lon'] for wf in wind_farm_list]
map_center_lat = sum(latitudes) / len(latitudes)
map_center_lon = sum(longitudes) / len(longitudes)

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.lines as mlines


country_style = {
    'UK': 'darkblue',
    'Norway': 'red',
    'Germany': 'green',
    'Denmark': 'blue',
    'Netherlands': 'orange'
}
# Create the Cartopy figure and axis.

fig_map = plt.figure(figsize=(14, 10))
ax_map = plt.axes(projection=ccrs.PlateCarree())
ax_map.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='black')
ax_map.add_feature(cfeature.OCEAN, facecolor='whitesmoke')
ax_map.set_extent([-5, 13, 50, 62], crs=ccrs.PlateCarree())

# A dictionary to collect one legend handle per country (to avoid duplicates).
legend_handles = {}

# Plot each wind farm with the color defined by its country.
for wf in wind_farm_list:
    country = wf['country']
    # Get the color; if country is not in the mapping, use gray.
    color = country_style.get(country, 'gray')
    
    # Plot the wind farm location.
    ax_map.plot(wf['lon'], wf['lat'], 'o', color=color, markersize=8,
                transform=ccrs.PlateCarree())
    
    # Create a legend handle if this country's marker isn't added yet.
    if country not in legend_handles:
        legend_handles[country] = mlines.Line2D([], [], color=color,
                                                marker='o', linestyle='None',
                                                markersize=8, label=country)

# Add the legend to the map.
ax_map.legend(handles=list(legend_handles.values()),
              loc='lower right', bbox_to_anchor=(0.99, 0.01), title="Wind Farms by Country")

# Optionally add a title and labels.
ax_map.set_title("Wind Farms Map by Country")
ax_map.set_xlabel("Longitude")
ax_map.set_ylabel("Latitude")

plt.savefig("wind_farms_map.png", dpi=150, bbox_inches="tight", pad_inches=0.1)

# Display the map.
plt.show()

# Define marker styles based on country
country_style = {
    'UK': {'color': 'darkblue', 'icon': 'bolt'},
    'Norway': {'color': 'red', 'icon': 'bolt'},
    'Germany': {'color': 'green', 'icon': 'bolt'},
    'Denmark': {'color': 'blue', 'icon': 'bolt'},
    'Netherlands': {'color': 'orange', 'icon': 'bolt'}
}



out_file = Path("wind_farms_final.txt")
cols      = ["name",
             "lat",
             "lon",
             "installed_capacity",
             "development_status",
             "country"]

# Write a tab-separated file
pruned_gdf[cols].to_csv(
    out_file,
    sep="\t",
    float_format="%.6f",
    index=False,
    lineterminator="\n"      
)

csv_file = Path("wind_farms_final.csv")
pruned_gdf[cols].to_csv(
    csv_file,
    index=False,       
    float_format="%.6f",
    lineterminator="\n",
)

cols = [
    "name",
    "lat",
    "lon",
    "installed_capacity",
    "development_status",
    "country",
]

xlsx_file = Path("wind_farms_final.xlsx")
pruned_gdf[cols].to_excel(xlsx_file, index=False)

print(f"Wrote {len(pruned_gdf)} rows to {xlsx_file.resolve()}")

print(f"Wrote {len(pruned_gdf)} wind farms to {out_file.resolve()}")
