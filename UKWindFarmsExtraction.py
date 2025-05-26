import pandas as pd
from pyproj import CRS, Transformer

# Load the CSV file with appropriate encoding and delimiter
file_path = '/Users/markusagersborg/Library/CloudStorage/OneDrive-NTNU/Prosjektoppgave/Optimeringmodell/UKOffshoreWindFarmsQ2_2024.csv'
wind_farms_df = pd.read_csv(file_path, encoding="ISO-8859-1", delimiter=";")

# Convert the 'Installed Capacity (MWelec)' column to a numeric type after handling commas
wind_farms_df['Installed Capacity (MWelec)'] = pd.to_numeric(
    wind_farms_df['Installed Capacity (MWelec)'].str.replace(',', '.'), errors='coerce'
)

# Filter for wind farms that have coordinate data (non-NaN X and Y coordinates)
wind_farms_with_coords = wind_farms_df.dropna(subset=['X-coordinate', 'Y-coordinate'])

# Set up coordinate systems for conversion: British National Grid to WGS84
bng_proj = CRS("EPSG:27700")  # British National Grid
wgs84_proj = CRS("EPSG:4326")  # WGS84
transformer = Transformer.from_crs(bng_proj, wgs84_proj, always_xy=True)

# Filtering criteria
min_capacity = 0  # Minimum installed capacity in MWelec
#Exclude wind farms on west side of UK
lat_condition = 55  
long_condition = -1.5 

# Convert coordinates and filter based on installed capacity, geographic bounds, and exclusion conditions
UK_filtered_wind_farms = []
for _, row in wind_farms_with_coords.iterrows():
    lon, lat = transformer.transform(row['X-coordinate'], row['Y-coordinate'])
    # Exclude wind farms that have both latitude < 55 and longitude < -1.5
    if (row['Installed Capacity (MWelec)'] >= min_capacity and
            not (lat < lat_condition and lon < long_condition)):
        UK_filtered_wind_farms.append({
            'name': row['Site Name'],
            'lat': round(lat, 5),
            'lon': round(lon, 5),
            'installed_capacity': row['Installed Capacity (MWelec)'],
            'development_status': row['Development Status'],
            'country': 'UK'
        })
