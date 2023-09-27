import pyproj
import matplotlib.pyplot as plt


# this class gives methods of transforming the longitude and latitude into coordination under CGCS2000
class CoordinationTransformer:
    def __init__(self, ref_longitude, ref_latitude):
        # Define the WGS84 coordinate reference system
        wgs84 = pyproj.CRS("EPSG:4326")

        # Define the Shanghai local CRS (e.g., EPSG:4549)
        # China Geodetic Coordinate System 2000 / 3-degree Gauss-Kruger CM 120E
        local_crs = pyproj.CRS("EPSG:4549")
        self.transformer_geo_to_cart = pyproj.Transformer.from_crs(
            wgs84,
            local_crs
        )
        self.transformer_cart_to_geo = pyproj.Transformer.from_crs(
            local_crs,
            wgs84
        )
        self.ref_longitude = ref_longitude
        self.ref_latitude = ref_latitude
        self.y, self.x = self.transformer_geo_to_cart.transform(ref_latitude, ref_longitude)
        # self.x, self.y = 0, 0

    def geoToCart(self, longitude, latitude):
        y, x = self.transformer_geo_to_cart.transform(latitude, longitude)
        return x - self.x, y - self.y

    def cartToGeo(self, x, y):
        lat, lon = self.transformer_cart_to_geo.transform(y + self.y, x + self.x)
        return lon, lat



if __name__ == "__main__":
    # 121.4294481, 31.0220344 is the lon, lat of the SJTU nan soccer field
    transformer = CoordinationTransformer(121.4294481, 31.0220344)
    with open('routes/trajectory.csv', 'r') as file:
        lines = file.readlines()

    x_coordinates = []
    y_coordinates = []

    # read file and record the lonlat inside it transform it into CGCS2000 coordination
    for line in lines:
        parts = line.strip().split(',')
        if len(parts) == 2:
            x, y = map(float, parts)
            x, y = transformer.geoToCart(x, y)
            x_coordinates.append(x)
            y_coordinates.append(y)

    print(((x_coordinates[2]-x_coordinates[1])**2+(y_coordinates[2]-y_coordinates[1])**2)**(1/2))

    plt.scatter(x_coordinates, y_coordinates, color='blue', label='Points')

    # title and tag
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Scatter Plot of Coordinates')

    # show figure
    plt.legend()
    plt.grid(True)
    plt.show()

    for i in range(3000):
        print(f'{x_coordinates[i]},{y_coordinates[i]}')


