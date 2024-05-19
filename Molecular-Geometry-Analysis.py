#MolMolecular-Geometry-Analysis
# 1 Bond Length 
# 2 Bond angle 
# 3 Out-of-Plane Angles
# 4 Torsion/Dihedral Angles
# 5 Center-of-Mass Translation
# 6 Principal Moments of Inertia
# 7 Rotational Constants
import math
import numpy as np 
an2masses = {
    1: 1.00784,
    2: 4.002602,
    3: 6.94,
    4: 9.0122,
    5: 10.81,
    6: 12.011,
    7: 14.007,
    8: 15.999,
    9: 18.998,
    10: 20.180,
}
def unit(geom, cart, a, b):
    if a == b:
        return 0.0
    else:
        return -(geom[a][cart + 1] - geom[b][cart + 1]) / bond_length(geom, a, b)
        
def bond_length(geom, a, b):
    if a == b:
        return 0.0
    else:
        return math.sqrt((geom[a][1] - geom[b][1]) ** 2 +
                         (geom[a][2] - geom[b][2]) ** 2 +
                         (geom[a][3] - geom[b][3]) ** 2)
        
def angle(geom, a, b, c):
    if a == b or b == c or a == c:
        return 0.0
    else:
        u0 = unit(geom, 0, b, a)
        u1 = unit(geom, 1, b, a)
        u2 = unit(geom, 2, b, a)
        v0 = unit(geom, 0, b, c)
        v1 = unit(geom, 1, b, c)
        v2 = unit(geom, 2, b, c)
        dot_product = u0 * v0 + u1 * v1 + u2 * v2
        if abs(dot_product) > 1.0:
            dot_product = 1.0 if dot_product > 0 else -1.0
        return math.acos(dot_product)

def outer_angle(geom, a, b, c, d):
    ebcd_x = (unit(geom, 1, c, b) * unit(geom, 2, c, d) - unit(geom, 2, c, b) * unit(geom, 1, c, d))
    ebcd_y = (unit(geom, 2, c, b) * unit(geom, 0, c, d) - unit(geom, 0, c, b) * unit(geom, 2, c, d))
    ebcd_z = (unit(geom, 0, c, b) * unit(geom, 1, c, d) - unit(geom, 1, c, b) * unit(geom, 0, c, d))
    exx = ebcd_x * unit(geom, 0, c, a)
    eyy = ebcd_y * unit(geom, 1, c, a)
    ezz = ebcd_z * unit(geom, 2, c, a)
    if math.sin(angle(geom, b, c, d)) != 0.0:
        theta = (exx + eyy + ezz) / math.sin(angle(geom, b, c, d))
    else:
        return 0.0
    if theta < -1.0:
        theta = math.asin(-1.0)
    elif theta > 1.0:
        theta = math.asin(1.0)
    else:
        theta = math.asin(theta)
    return theta

def torsion_angle(geom, a, b, c, d):
    eabc_x = (unit(geom, 1, b, a) * unit(geom, 2, b, c) - unit(geom, 2, b, a) * unit(geom, 1, b, c))
    eabc_y = (unit(geom, 2, b, a) * unit(geom, 0, b, c) - unit(geom, 0, b, a) * unit(geom, 2, b, c))
    eabc_z = (unit(geom, 0, b, a) * unit(geom, 1, b, c) - unit(geom, 1, b, a) * unit(geom, 0, b, c))
    ebcd_x = (unit(geom, 1, c, b) * unit(geom, 2, c, d) - unit(geom, 2, c, b) * unit(geom, 1, c, d))
    ebcd_y = (unit(geom, 2, c, b) * unit(geom, 0, c, d) - unit(geom, 0, c, b) * unit(geom, 2, c, d))
    ebcd_z = (unit(geom, 0, c, b) * unit(geom, 1, c, d) - unit(geom, 1, c, b) * unit(geom, 0, c, d))
    exx = eabc_x * ebcd_x
    eyy = eabc_y * ebcd_y
    ezz = eabc_z * ebcd_z
    if( math.sin(angle(geom, a, b, c)) * math.sin(angle(geom, b, c, d)) != 0.0):
        tau = (exx + eyy + ezz) / (math.sin(angle(geom, a, b, c)) * math.sin(angle(geom, b, c, d)))
    else:
        return 0.0
    if tau < -1.0:
        tau = math.acos(-1.0)
    elif tau > 1.0:
        tau = math.acos(1.0)
    else:
        tau = math.acos(tau)

    cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y
    cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z
    cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x
    norm = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z
    if (norm!=0.0):
        cross_x /= norm
        cross_y /= norm
        cross_z /= norm
    sign = 1.0
    dot = cross_x * unit(geom, 0, b, c) + cross_y * unit(geom, 1, b, c) + cross_z * unit(geom, 2, b, c)
    if dot < 0.0:
        sign = -1.0

    return tau * sign

with open('data1.dat', 'r') as input_file:
    natom = int(input_file.readline().strip())
    geom = []
    for i in range(natom):
        line = input_file.readline().strip()
        data = line.split()
        zval = int(data[0])
        x = float(data[1])
        y = float(data[2])
        z = float(data[3])
        geom.append((zval, x, y, z))
    print(f"Number of atoms: {natom}")
    print("Input Cartesian coordinates:")
    for atom in geom:
        print(f"({atom[0]}) ({atom[1]:20.12f}) ({atom[2]:20.12f}) ({atom[3]:20.12f})")

    for j in range(natom):
        for k in range(j):
            R = bond_length(geom, j, k)
            print(j, k, R)

    for i in range(natom):
        for j in range(natom):
            for k in range(natom):
                ang_degrees = angle(geom, i, j, k) * (180.0 / math.pi)
                print(i, j, k, ang_degrees)

    for i in range(natom):
        for j in range(natom):
            for k in range(natom):
                for l in range(j):
                    dihedral_degrees = outer_angle(geom, i, j, k, l) * (180.0 / math.pi)
                    print(i, j, k, l, dihedral_degrees)

    for i in range(natom):
        for j in range(i):
            for k in range(j):
                for l in range(k):
                    torsional_angle = torsion_angle(geom, i, j, k, l) * (180 / math.pi)
                    print(i, j, k, l, torsional_angle)

    M = 0.0
    for i in range(natom):
        M += an2masses[int(geom[i][0])]  # Corrected line

    xcm = 0.0
    ycm = 0.0
    zcm = 0.0
    for i in range(natom):
        mi = an2masses[int(geom[i][0])] 
        xcm += mi * geom[i][1]
        ycm += mi * geom[i][2]
        zcm += mi * geom[i][3]
    xcm /= M
    ycm /= M
    zcm /= M

    print(f"Xcm: {xcm}")
    print(f"Ycm: {ycm}")
    print(f"Zcm: {zcm}")
    I = np.zeros((3,3))
    for i in range(natom):
        mi = an2masses[int(geom[i][0])]
        x = geom[i][1] - xcm
        y = geom[i][2] - ycm
        z = geom[i][3] - zcm
        I[0][0] += mi * (y**2 + z**2)
        I[1][1] += mi * (x**2 + z**2)
        I[2][2] += mi * (x**2 + y**2)
        I[0][1] -= mi * x * y
        I[0][2] -= mi * x * z
        I[1][2] -= mi * y * z
    I[1][0] = I[0][1]
    I[2][0] = I[0][2]
    I[2][1] = I[1][2]

    print(f"Moment of inertia tensor :\n{I}")

    evals, evecs = np.linalg.eig(np.array(I))
    evals = np.sort(evals)  

    print(f"Principal moments of inertia:\n{evals}")

    amu_bohr2_to_amu_AA2 = 0.529177249**2
    principal_moments_AA2 = evals * amu_bohr2_to_amu_AA2
    print(f"Principal moments of inertia (amu AA^2)\n{principal_moments_AA2}")
    
    amu_bohr2_to_g_cm2 = 1.6605402E-24 * 0.529177249E-8**2
    principal_moments_g_cm2 = evals * amu_bohr2_to_g_cm2
    print(f"Principal moments of inertia (g cm^2)\n{principal_moments_g_cm2}")

    if natom == 2:
        print("Molecule is diatomic")
    elif evals[0] < 1e-4:
        print("Molecule is linear")
    elif np.isclose(evals[0], evals[1]) and np.isclose(evals[1], evals[2]):
        print("Molecule is a spherical top")
    elif np.isclose(evals[0], evals[1]):
        print("Molecule is an oblate symmetric top")
    elif np.isclose(evals[1], evals[2]):
        print("Molecule is a prolate symmetric top")
    else:
        print("Molecule is an asymmetric top")


import numpy as np
h = 6.6260755E-34  
m = 1.6605402E-27  
r = 0.529177249E-10 
c = 2.99792458E10  
I = m * (r ** 2)  
B1 = h / (8.0 * np.pi * np.pi * I * c)   
A = B1 / evals[0]
B = B1 / evals[1]
C = B1 / evals[2]
print("\nRotational constants (cm⁻¹):")
print(f"\tA = {A:.6e} cm⁻¹\t B = {B:.6e} cm⁻¹\t C = {C:.6e} cm⁻¹")

