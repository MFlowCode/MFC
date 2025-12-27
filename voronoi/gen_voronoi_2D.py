import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import freud


# lloyd relaxation
def compute_simplex_centroid(simplex_vertices):
    v1 = simplex_vertices[:, :, 0]
    v2 = simplex_vertices[:, :, 1]
    v3 = simplex_vertices[:, :, 2]

    v1_mean = np.mean(v1, axis=1)
    v2_mean = np.mean(v2, axis=1)
    v3_mean = np.mean(v3, axis=1)

    simplex_centroids = np.array([v1_mean, v2_mean, v3_mean])

    return simplex_centroids

def compute_simplex_area(simplex_vertices):
    v1 = simplex_vertices[:, :, 0]
    v2 = simplex_vertices[:, :, 1]
    v3 = simplex_vertices[:, :, 2]

    area = 0.5 * np.linalg.norm( np.cross(v2 - v1, v3 - v1), axis=1 )

    return area

def lloyd_relaxation_2d(initial_points, box, w=1.0, iterations=20):
    points = initial_points

    for _ in range(iterations):
        voro = freud.locality.Voronoi()
        voro_data = voro.compute((box, initial_points))
        vertices = voro_data.polytopes

        for i in range(len(points)):
            n = len(vertices[i])

            simplex_vertices = np.array( [(points[i, :], vertices[i][j-1], vertices[i][j]) for j in range(n)] )

            simplex_centroids = compute_simplex_centroid(simplex_vertices)
            simplex_areas = compute_simplex_area(simplex_vertices)

            centroid = (1/np.sum(simplex_areas)) * (np.sum(simplex_centroids*simplex_areas, axis=1))

            dist = centroid - points[i, :]

            points[i, :] += w * dist

        points = box.wrap(points)

    return points

if (__name__ == '__main__'):
    print('running 2D...')

    # setup 
    phi = 0.4
    D = 0.1
    L = 10*D

    N = int( 4*phi*L**2 / (np.pi*D**2) )
    print(f'volume fraction phi: {phi}, number of circles: {N}')

    x_i = L/2 * np.random.uniform(-1, 1, N)
    y_i = L/2 * np.random.uniform(-1, 1, N)
    z_i = L/2 * np.random.uniform(-1, 1, N) * 0

    initial_points = np.stack((x_i, y_i, z_i), axis=1)

    box = freud.box.Box.square(L)
    voro = freud.locality.Voronoi()

    cells = voro.compute((box, initial_points)).polytopes

    # plot initial distribution
    plt.figure()
    ax = plt.gca()
    voro.plot(ax=ax, cmap='RdBu')
    ax.scatter(initial_points[:, 0], initial_points[:, 1], s=5, c='k')
    plt.show()
    plt.close()

    # calculate relaxed points
    relaxed_points = lloyd_relaxation_2d(initial_points, box, w=1.5, iterations=25)
    voro.compute((box, relaxed_points))

    # plot relaxed distribution
    plt.figure()
    ax = plt.gca()
    voro.plot(ax=ax, cmap='RdBu')
    ax.scatter(relaxed_points[:, 0], relaxed_points[:, 1], s=5, c='k')
    plt.show()
    plt.close()



