"""
Calculate the average annual storm coverage ratio.
Ratio = (Average total storm area per year) / (Total region area)
"""
import json
from pathlib import Path
from shapely.geometry import shape, box
from shapely.ops import transform
import pyproj

# Load storm data
storms_by_year = {}
hail_maps_dir = Path('hail_maps')

for year in range(2011, 2021):
    with open(hail_maps_dir / f'hail_{year}.geojson', 'r') as f:
        data = json.load(f)
    storms_by_year[year] = [shape(feature['geometry']) for feature in data['features']]

# Get region bounds from all storms
min_lon = min_lat = float('inf')
max_lon = max_lat = float('-inf')

for storms in storms_by_year.values():
    for storm in storms:
        bounds = storm.bounds
        min_lon = min(min_lon, bounds[0])
        min_lat = min(min_lat, bounds[1])
        max_lon = max(max_lon, bounds[2])
        max_lat = max(max_lat, bounds[3])

# Set up projection for accurate area calculation
wgs84 = pyproj.CRS('EPSG:4326')
albers = pyproj.CRS('ESRI:102003')
#project is a function that takes in x,y in lon, lat and returns x,y in meters
project = pyproj.Transformer.from_crs(wgs84, albers, always_xy=True).transform

# Calculate region area
# print(min_lon, min_lat, max_lon, max_lat)
region_box = box(min_lon, min_lat, max_lon, max_lat)
region_area = transform(project, region_box).area # have to convert before square degrees is meaningless

# Calculate total storm area for each year
yearly_storm_areas = []
for year, storms in storms_by_year.items():
    total_area = sum(transform(project, storm).area for storm in storms)
    yearly_storm_areas.append(total_area)
    print(f"{year}: {total_area / 1e6:.2f} km²")

# Calculate average and ratio
avg_storm_area = sum(yearly_storm_areas) / len(yearly_storm_areas)
ratio = avg_storm_area / region_area

# Print results
print(f"\n{'='*60}")
print(f"Region bounds: ({min_lon:.2f}, {min_lat:.2f}, {max_lon:.2f}, {max_lat:.2f})")
print(f"Region area: {region_area / 1e6:.2f} km²")
print(f"Average annual storm area: {avg_storm_area / 1e6:.2f} km²")
print(f"Coverage ratio: {ratio:.6f} ({ratio*100:.4f}%)")
print(f"{'='*60}")