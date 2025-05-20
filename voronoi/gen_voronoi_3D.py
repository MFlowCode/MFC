import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import freud


# lloyd relaxation
def compute_tetrahedron_centroid(tetrahedron_vertices):
    
    return np.mean(tetrahedron_vertices, axis=0)

def compute_tetrahedron_volume(tetrahedron_vertices):
    v0, v1, v2, v3 = tetrahedron_vertices
    matrix = np.vstack([v1 - v0, v2 - v0, v3 - v0]).T
    volume = np.abs(np.linalg.det(matrix)) / 6

    return volume

def lloyd_relaxation_3d(initial_points, box, w=1, iterations=10):
    points = initial_points

    for _ in range(iterations):
        voro = freud.locality.Voronoi()
        voro_data = voro.compute((box, points))
        vertices = voro_data.polytopes

        for i in range(len(points)):
            n = len(vertices[i])

            tetrahedra = []
            for j in range(n):
                tetrahedra.append([points[i, :], vertices[i][j], vertices[i][(j+1) % n], vertices[i][(j+2) % n]])

            centroids = np.array([compute_tetrahedron_centroid(t) for t in tetrahedra])
            volumes = np.array([compute_tetrahedron_volume(t) for t in tetrahedra])

            weighted_centroid = np.sum(centroids * volumes[:, np.newaxis], axis=0)
            total_volume = np.sum(volumes)

            if total_volume > 1.0e-12:
                centroid = weighted_centroid / total_volume
                dist = centroid - points[i, :]

                points[i, :] += w * dist

        points = box.wrap(points)

    return points

if (__name__ == '__main__'): 
    print('running 3D...')

    # setup 
    phi = 0.05
    str_phi = '005'

    D = 0.1
    L = 10*D

    output_dir = '../examples/phi'+str_phi
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    N_sphere = int( 6*phi*L**3 / (np.pi*D**3) )
    print(f'volume fraction phi: {phi}, number of spheres: {N_sphere}')
    print(f'actual phi value: {N_sphere*4/3*np.pi*(D/2)**3/(L**3)}')

    x_i = L/2 * np.random.uniform(-1, 1, N_sphere)
    y_i = L/2 * np.random.uniform(-1, 1, N_sphere)
    z_i = L/2 * np.random.uniform(-1, 1, N_sphere)

    initial_points = np.stack((x_i, y_i, z_i), axis=1)
    box = freud.box.Box.cube(L)
    
    relaxed_points = lloyd_relaxation_3d(initial_points, box, iterations=30)
    print(np.shape(relaxed_points))

    np.savetxt(output_dir+'/sphere_array_locations.txt', relaxed_points)

    # check no spheres are overlaping
    for i in range(N_sphere):
        for j in range(N_sphere):
            if (i != j):
                dist = np.sqrt((relaxed_points[i, 0] - relaxed_points[j, 0])**2 + (relaxed_points[i, 1] - relaxed_points[j, 1])**2 + (relaxed_points[i, 2] - relaxed_points[j, 2])**2)
                if (dist <= 1.05*D):
                    print(f'spheres overlaping, dist={dist}, spheres #: {i}, {j}')
                    print(f'locations: ({relaxed_points[i, :]}), ({relaxed_points[j, :]})')

    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(initial_points[:, 0], initial_points[:, 1], initial_points[:, 2], color='blue', s=10)
    ax1.set_title('initial points')
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(relaxed_points[:, 0], relaxed_points[:, 1], relaxed_points[:, 2], color='red', s=10)
    ax2.set_title('relaxed points')
    plt.show()
    plt.close()
