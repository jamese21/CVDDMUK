from shapely.geometry import Polygon, MultiPolygon, Point
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import random
# import folium

# method for plotting a polygon or multipolygon
def plot_map(geometry):
    # Ensure the input is a list for consistency
    if isinstance(geometry, (Polygon, MultiPolygon)):
        geometry = [geometry]  # Wrap single geometries in a list

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(geometry=geometry)

    # Plotting
    fig, ax = plt.subplots()
    gdf.plot(ax=ax, alpha=0.5, edgecolor='black')

    # Set axis properties
    ax.set_title('Polygons and Multipolygons')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.axis('equal')  # Equal scaling
    plt.show()

# plot map and label each region in order
def plot_map_with_labels(geometry):
    # Ensure the input is a list for consistency
    if isinstance(geometry, (Polygon, MultiPolygon)):
        geometry = [geometry]  # Wrap single geometries in a list

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(geometry=geometry)

    # Plotting
    fig, ax = plt.subplots()
    gdf.plot(ax=ax, alpha=0.5, edgecolor='black')

    # Label each polygon with its order number
    for idx, geom in enumerate(gdf.geometry):
        # Calculate the centroid of each polygon
        centroid = geom.centroid
        # Plot the label at the centroid
        ax.text(centroid.x, centroid.y, str(idx + 1), fontsize=12, ha='center', va='center', color='red')

    # Set axis properties
    ax.set_title('Polygons and Multipolygons')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.axis('equal')  # Equal scaling
    plt.show()

# plot the area inside of a polygon
def plot_area(polygon):
    x, y = polygon.exterior.xy

    # Plot the polygon
    plt.figure()
    plt.fill(x, y, alpha=0.5, fc='blue', ec='black')  # Fill with blue color
    plt.title('Polygon from WKT')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.axis('equal')  # Keep aspect ratio
    plt.show()

# Keep only the largest polygon of each multipolygon (remove islands)
def remove_islands(geometry):
    no_islands = []
    for shape in geometry:
        # plot_map(shape)
        if shape.geom_type=='MultiPolygon':
            polygons = list(shape.geoms)
            max_size = 0
            max_poly = None
            for p in polygons:
                if p.area > max_size:
                    max_size = p.area
                    max_poly = p
            no_islands.append(max_poly)
        else:
            no_islands.append(shape)
    return no_islands

# Functions for generating points across polygon with different distribution
def generate_points_evenly(min_x, max_x, min_y, max_y, polygon):
    x_diff = max_x-min_x
    y_diff = max_y-min_y
    x_ratio = x_diff/(x_diff+y_diff)
    y_ratio = y_diff/(x_diff+y_diff)
    x_num = int(x_ratio*100) # could be missing edge
    y_num = int(y_ratio*100)
    x_gap = x_diff/x_num
    y_gap = y_diff/y_num

    points = []

    for x_count in range(x_num):
        for y_count in range(y_num):
            x_coord = min_x+x_gap*x_count
            y_coord = min_y+y_gap*y_count
            point = Point(x_coord, y_coord)
            if polygon.contains(point):
                points.append(point)
    return points

def generate_points_no_outline(min_x, max_x, min_y, max_y, polygon):
    density = 100
    x_diff = max_x-min_x
    y_diff = max_y-min_y
    x_ratio = x_diff/(x_diff+y_diff)
    y_ratio = y_diff/(x_diff+y_diff)
    x_num = int(x_ratio*density) # could be missing edge
    y_num = int(y_ratio*density)
    x_gap = x_diff/x_num
    y_gap = y_diff/y_num

    points = []

    for x_count in range(x_num):
        for y_count in range(y_num):
            x_coord = min_x+x_gap*x_count
            y_coord = min_y+y_gap*y_count
            point = Point(x_coord, y_coord)
            if polygon.contains(point):
                if point.distance(polygon.boundary) > density/5000:
                    points.append(point)
    return points

# Functions for generating points across polygon with different distribution
def place_n_points(polygon, n, iterations=5):
    # values
    if polygon == None:
        return []
    area = polygon.area
    r = (area / (n * np.pi)) ** 0.5

    # n random points poisson disc distribution (ish)
    min_x, min_y, max_x, max_y = polygon.bounds
    points = []
    while len(points) < n:
        x = random.uniform(min_x, max_x)
        y = random.uniform(min_y, max_y)
        candidate = Point(x, y)

        # Append candidate if within the polygon and distance from other points >= r
        if polygon.contains(candidate):
            if all(candidate.distance(point) >= r for point in points):
                points.append(candidate)
    return points

def distribute_points(polygon, points, iterations):
    if polygon == None:
        return []
    area = polygon.area
    r = (area / (len(points) * np.pi)) ** 0.5
    for _ in range(iterations):
        centroids = []
        for point in points:
            neighborhood = [p for p in points if p.distance(point) < 2*r]
            if neighborhood:
                centroid_x = sum(p.x for p in neighborhood) / len(neighborhood)
                centroid_y = sum(p.y for p in neighborhood) / len(neighborhood)
                centroids.append(Point(centroid_x, centroid_y))
        points = centroids
    return points

# plot a geometry and a set of points over it
def plot_points_and_map(points, geometry, point_size):
    # Ensure the input is a list for consistency
    if isinstance(geometry, (Polygon, MultiPolygon)):
        geometry = [geometry]  # Wrap single geometries in a list

    gdf_polygons = gpd.GeoDataFrame(geometry=geometry)
    gdf_points = gpd.GeoDataFrame(geometry=points)

    # Plotting
    fig, ax = plt.subplots()
    gdf_polygons.plot(ax=ax, alpha=0.5, edgecolor='black')
    gdf_points.plot(ax=ax, color='red', marker='.', markersize=point_size, label='Points')

    # Set axis properties
    ax.set_title('Polygons and Points')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.axis('equal')  # Equal scaling
    plt.show()

def plot_points_and_map_interactive(points, geometry):
    # Initialize a map centered around the first point (or a central location)
    if points:
        center = [points[0].y, points[0].x]
    else:
        center = [geometry[0].centroid.y, geometry[0].centroid.x]
    m = folium.Map(location=center, zoom_start=12)

    # Add geometries (polygons)
    for geom in geometry:
        folium.GeoJson(data=mapping(geom), style_function=lambda x: {
            'fillColor': 'blue',
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.5
        }).add_to(m)

    # Add points
    for point in points:
        folium.CircleMarker(location=[point.y, point.x], radius=3, color='red').add_to(m)

    return m

# Convert a list of polygons/multipolygons to a single multipolygon
def list_to_multipolygon(geom_list):
    all_polygons = []
    for geom in geom_list:
        if isinstance(geom, Polygon):
            all_polygons.append(geom)
        elif isinstance(geom, MultiPolygon):
            all_polygons.extend(geom.geoms)
    merged = MultiPolygon(all_polygons)
    return merged