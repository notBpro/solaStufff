import math
import json
import numpy as np
from pathlib import Path
from shapely.geometry import Point, shape
from shapely.prepared import prep
from shapely.ops import transform
import pyproj

# Min distance ensures d=0 for points inside storm polygons (actual impact),
# while blurring uncertainty for points outside (potential future impact).

def min_distance_m(point: Point, poly) -> float:
    """
    Calculate the minimum distance in meters from a Shapely Point to a polygon.
    
    point: shapely.geometry.Point
        The point to measure from (must be in a projected CRS with meter units)
    poly: shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        The polygon to measure to (must be in same CRS as point)
    
    Returns:
        Distance in meters. Returns 0 if the point is inside the polygon.
    """
    # Check if point is inside the polygon
    prepared_poly = prep(poly) #makes contains log n instead checking all edges
    if prepared_poly.contains(point):
        return 0.0
    
    # Calculate distance (already in meters due to projected CRS)
    return point.distance(poly)

# def centroid_distance_m(point: Point, poly) -> float:
#     """
#     Calculate the distance in meters from a point to the centroid of a polygon.
    
#     Args:
#         point: shapely.geometry.Point (must be in projected CRS with meter units)
#         poly: shapely.geometry.Polygon or MultiPolygon (must be in same CRS as point)
    
#     Returns:
#         Distance in meters from point to polygon's centroid
#     """
#     centroid = poly.centroid
#     return point.distance(centroid)

def calculate_storm_probability(storms_by_year, p, h):
    """
    Calculate storm probability at point p using Gaussian kernel smoothing.
    
    Args:
        storms_by_year: dict mapping year -> list of storm polygons with their areas
        [year: [storm_poly1, storm_poly2, ...]]
        p: tuple (latitude, longitude) 
        h: blur distance in meters (the knob)
    
    Returns:
        P_final: final probability averaged across all years
    """
    
    # Kernel function: normalized Gaussian - dome shape
    def kernel(d, h):
        if h == 0:
            return 1.0 if d == 0 else 0.0
        return math.exp(-(d**2)/(2*h**2)) / (2 * math.pi * h**2) # e^(-d^2 / 2h^2) / (2 * pi * h^2)

    
    total_prob = 0
    
    # Get number of years for averaging
    number_of_years = len(storms_by_year)
    
    # For each year Y in storms_by_year
    for year, storms in storms_by_year.items():
        # Initialize yearly sum (total "paint" from all storms this year)
        yearly_sum = 0

        # For each storm S in year Y
        for storm_poly in storms:
            # Get area of storm S in m²
            # Use Albers Equal Area projection for accurate area calculation
            wgs84 = pyproj.CRS('EPSG:4326')
            albers = pyproj.CRS('ESRI:102003') 

            #project is a function that takes in x,y in lon, lat and returns x,y in meters
            project = pyproj.Transformer.from_crs(wgs84, albers, always_xy=True).transform
            
            # Project polygon to equal area and calculate area in square meters
            storm_projected = transform(project, storm_poly) #storm_poly looks like this [[[x1,y1], [x2,y2], ...], [[x1,y1], [x2,y2], ...], ...]
            area_S = storm_projected.area
            
            # Convert point p to Shapely Point
            point = Point(p[0], p[1])  # (lon, lat) for Shapely
            
            # Project both point and polygon to meter-based coordinate system for distance calculation
            point_projected = transform(project, point)
            storm_projected = transform(project, storm_poly)
            
            # Calculate distance between p and S (0 if inside) - now in meters
            d = min_distance_m(point_projected, storm_projected)
            
            # Calculate influence: how much of storm S's paint lands on point p
            influence = area_S * kernel(d, h)
            
            # Add to yearly sum
            yearly_sum += influence
        
        # Convert yearly_sum to probability (cap at 1 only 1 payout per year)
        # cannot be more certain than certain
        P_Y = min(1, yearly_sum) #prob for year y

        # Add to total probability
        total_prob += P_Y
    
    # Average across all years
    P_final = total_prob / number_of_years
    
    return P_final

if __name__ == "__main__":
    # Load storm data from your hail maps
    storms_by_year = {}
    hail_maps_dir = Path('hail_maps')
    
    for year in range(2011, 2021):
        with open(hail_maps_dir / f'hail_{year}.geojson', 'r') as f: #combines the directory with the filename to open it and read
            data = json.load(f) #json to python object
        
        # Extract storm polygons from GeoJSON and store in storms_by_year
        storms = [shape(feature['geometry']) for feature in data['features']]
        storms_by_year[year] = storms
    
    # # Define region bounds (region data taken from baseline)
    min_lat, max_lat = 25.77, 36.66
    min_lon, max_lon = -106.61, -93.35

    num_points = 25
    
    # Generate random points within the region(uniform for every point gets equal prob of being picked)
    np.random.seed(69)  # For reproducibility
    random_points = [(np.random.uniform(min_lon, max_lon), 
                      np.random.uniform(min_lat, max_lat)) 
                     for _ in range(num_points)]
    
    # # Blur distance in meters
    h = 30000  # 30km blur, larger h means more blurring
    
    # Calculate probability for each point and average
    probabilities = []
    for i, p in enumerate(random_points):
        probability = calculate_storm_probability(storms_by_year, p, h)
        probabilities.append(probability)
        print(f"Point {i+1}/{num_points} at {p}: {probability:.6f}")
    
    # Calculate and print average
    avg_probability = np.mean(probabilities)
    print(f"\nAverage probability across {num_points} random points: {avg_probability:.6f}")

    point = (-96.7838, 32.7839)
    probability = calculate_storm_probability(storms_by_year, point, h)
    print(f"Probability at point {point}: {probability:.6f}")