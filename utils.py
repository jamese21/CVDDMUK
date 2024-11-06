from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd
import matplotlib.pyplot as plt

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

# plot a geometry and a set of points over it
def plot_points_and_map(points, geometry):
    # Ensure the input is a list for consistency
    if isinstance(geometry, (Polygon, MultiPolygon)):
        geometry = [geometry]  # Wrap single geometries in a list

    gdf_polygons = gpd.GeoDataFrame(geometry=geometry)
    gdf_points = gpd.GeoDataFrame(geometry=points)

    # Plotting
    fig, ax = plt.subplots()
    gdf_polygons.plot(ax=ax, alpha=0.5, edgecolor='black')
    gdf_points.plot(ax=ax, color='red', marker='.', markersize=10, label='Points')

    # Set axis properties
    ax.set_title('Polygons and Points')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.axis('equal')  # Equal scaling
    plt.show()