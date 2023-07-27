    # def find_backbone_orientation(self):
    #     """
    #     This function returns the relative orientation of two helicenes (aromatic molecules) relative to each other, based on the 
    #     aromatic rings, to check 
    #     whether the backbones are parallel, antiparallel, or orthogonal to each other. 

    #     The orientation is computed based on the centroid and the midpoint between the two terminal rings of each molecule.
    #     """
    #     rings_per_mol = int(int(len(self.rings))/2)

    #     mola = np.array(self.centres[:rings_per_mol])
    #     molb = np.array(self.centres[rings_per_mol:])
    #     assert euclidean3d(self.centreA, mola.mean(axis=0)) < euclidean3d(self.centreB, mola.mean(axis=0)) 
    #     # define molecule vector orthogonal to disk plane along which to measure the orientation.

# orientation functions/tools/helpers:
# def check_dist_from_centroid(atom, centroid):
#     dist = euclidean3d(atom, centroid)
#     return np.round(dist, 2)

# def midpoint(p1, p2):
#     x1, y1, z1 = p1
#     x2, y2, z2 = p2
#     x_mid = (x1 + x2) // 2
#     y_mid = (y1 + y2) // 2
#     z_mid = (z1 + z2) // 2
#     return [x_mid, y_mid, z_mid]