# 3d_Curved_Beam_Element_with_varying_reinforcement_and _stirups

import openseespy.opensees as ops
import math
import matplotlib.pyplot as plt

# Define material properties for concrete and steel
Ec = 26600.0  # Young's modulus of concrete in MPa
nc = 0.2      # Poisson's ratio for concrete
rc = 2.5e-9   # Density of concrete in tonnes/mm^3 (2500 kg/m^3)
fc = 30.0     # Compressive strength of concrete in MPa
cs = -0.2 * fc  # Crushing strength of concrete

# Properties of grade 60 steel in MPa
Es = 200000.0  # Young's modulus for steel in MPa
fs = 413.60    # Yield strength of steel in MPa

# Define beam section dimensions
bw = 250.0     # Beam width in mm
bh = 500.0     # Beam height in mm
cv = 50.0      # Cover to reinforcement in mm

# Torsional properties using an approximate formula for the rectangular beam
b = bw
h = bh
J = (b * h**3 / 3) * ((16 / 3) - (3.36 * b / h) * (1 - (b**4 / (12 * h**4))))

# Shear modulus for concrete
Gc = Ec / (2 * (1 + nc))  # Shear modulus in MPa

# Prompt user for geometry parameters
r = float(input("Enter the radius of the circle (in mm, e.g., 1000): "))  # Radius of the circle
cp = input("Enter the circle portion like ('1', '1/2', '1/4'): ").strip().lower()  # Circle portion

# Validate circle portion input and determine angle range
if cp == '1/4':
    ar = math.pi / 2
elif cp == '1/2':
    ar = math.pi
elif cp == '1':
    ar = 2 * math.pi
else:
    raise ValueError("Invalid circle portion. Choose '1', '1/2', or '1/4'.")

# Calculate length of the curved beam
bl = r * ar
print(f"Calculated beam length: {bl:.2f} mm")

# Prompt user for uniformly distributed load values
ul = float(input("Enter the uniformly distributed load (Tonn/mm, e.g., 0.01): "))  # Uniform load

# Prompt user for number of integration points
np = int(input("Enter the number of integration points for the curve: "))  # Number of integration points

# Define reinforcement positions and diameters based on user input
rb_pos = []  # List to store rebar positions
rb_areas = []  # List to store rebar areas
nr = int(input("Enter the number of longitudnal rebars to be placed in the beam: "))  # Number of rebars
for _ in range(nr):
    coord_input = input(f"Enter the y,z coordinates for rebar {_+1} (origin of coordinate is crossection center, format: y,z): ")  # Rebar coordinates
    dia_input = float(input(f"Enter the diameter of rebar {_+1} (in mm): "))  # Rebar diameter
    y_pos, z_pos = map(float, coord_input.split(','))  # Parse coordinates
    rb_pos.append((y_pos, z_pos))  # Append positions
    rb_area = math.pi * (dia_input / 2)**2  # Calculate rebar area
    rb_areas.append(rb_area)  # Append areas

# Prompt user for stirrup properties
st_pos = []  # List to store stirrup positions
st_areas = []  # List to store stirrup areas
cl = 0.0  # Current length
while cl < bl:
    sp = float(input(f"Enter stirrup spacing from one edge(mm) (current total length: {cl:.2f} mm, max: {bl:.2f} mm): "))  # Stirrup spacing
    if cl + sp > bl:
        print("Stirrups provided based on your values. Proceeding to analysis.")
        break
    else:
        st_pos.append(cl)  # Append stirrup positions
        cl += sp  # Update current length
    dia_input = float(input(f"Enter the diameter of stirrup (in mm): "))  # Stirrup diameter
    st_area = math.pi * (dia_input / 2)**2  # Calculate stirrup area
    st_areas.append(st_area)  # Append areas

print(f"Total stirrup spacing: {sum(st_pos):.2f} mm")

# Initialize the OpenSees model
ops.wipe()  # Clear previous model
ops.model('basic', '-ndm', 3, '-ndf', 6)  # Define 3D model with 6 DOF per node
ops.geomTransf('Linear', 1, 0, 0, 1)  # Define geometric transformation

# Material properties of each fiber assigned as uniaxial stress-strain relation
ops.uniaxialMaterial('Concrete01', 1, -fc, -0.002, cs, -0.006)  # Concrete material
ops.uniaxialMaterial('Steel01', 2, fs, Es, 0.01)  # Steel material
ops.uniaxialMaterial('Elastic', 3, Gc)  # Torsional material
ops.uniaxialMaterial('Steel01', 4, fs, Es, 0.01)  # Stirrup material

# Define the section with fibers, variable rebar positions, and stirrups
ops.section('Fiber', 1, '-GJ', J)
# Patch for concrete core
ops.patch('rect', 1, 10, 10, -bw / 2 + cv, -bh / 2 + cv, bw / 2 - cv, bh / 2 - cv)
# Patch for concrete cover
ops.patch('rect', 1, 10, 1, -bw / 2, -bh / 2, bw / 2, -bh / 2 + cv)
ops.patch('rect', 1, 10, 1, -bw / 2, bh / 2 - cv, bw / 2, bh / 2)
ops.patch('rect', 1, 1, 10, -bw / 2, -bh / 2 + cv, -bw / 2 + cv, bh / 2 - cv)
ops.patch('rect', 1, 1, 10, bw / 2 - cv, -bh / 2 + cv, bw / 2, bh / 2 - cv)

# Adding longitudinal rebars
for i, pos in enumerate(rb_pos):
    ops.layer('straight', 2, 1, rb_areas[i], pos[0], pos[1], pos[0], pos[1])

# Adding stirrups
for i, pos in enumerate(st_pos):
    y1 = -bh / 2 + cv
    y2 = bh / 2 - cv
    z1 = -bw / 2 + cv
    z2 = bw / 2 - cv
    # Bottom stirrup - runs horizontally across the width at the bottom
    ops.layer('straight', 4, 2, st_areas[i], y1, z1, y1, z2)
    # Top stirrup - runs horizontally across the width at the top
    ops.layer('straight', 4, 2, st_areas[i], y2, z1, y2, z2)
    # Side stirrup 1 - runs vertically along one side of the beam
    ops.layer('straight', 4, 2, st_areas[i], y1, z1, y2, z1)
    # Side stirrup 2 - runs vertically along the other side of the beam
    ops.layer('straight', 4, 2, st_areas[i], y1, z2, y2, z2)

# Define nodes for the circle portion
for i in range(np + 1):
    th = ar * (i / np)
    x = r * math.cos(th)
    y = 0.0  # Y-coordinate set to 0.0 assuming the single segment
    z = r * math.sin(th)  # Z-coordinate (assuming the circle lies in the XZ-plane)
    ops.node(i, x, y, z)

# Fix the first node 0 to all 6 degrees of freedom in 3D study
ops.fix(0, 1, 1, 1, 1, 1, 1)

# Define 3D beam elements with fiber section between nodes
ops.beamIntegration('Lobatto', 1, 1, 3)
for i in range(np):
    ops.element('forceBeamColumn', i + 1, i, i + 1, 1, 1)

# Apply uniform load to the elements
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
for i in range(1, np + 1):
    # Length of each element is calculated here using the Euclidean approach
    l = math.sqrt((ops.nodeCoord(i)[0] - ops.nodeCoord(i-1)[0])**2 +
                  (ops.nodeCoord(i)[1] - ops.nodeCoord(i-1)[1])**2 +
                  (ops.nodeCoord(i)[2] - ops.nodeCoord(i-1)[2])**2)
    # Apply uniform load to the element
    ops.eleLoad('-ele', i, '-type', 'beamUniform', 0.0, 0.0, -ul * l, 0.0, 0.0, 0.0)

# Analysis techniques
ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
ops.test('EnergyIncr', 1e-5, 10000)
ops.algorithm('KrylovNewton')
ops.integrator('LoadControl', 0.00001)
ops.analysis('Static')

# Run the analysis incrementally
ns = 10000
for step in range(ns):
    if ops.analyze(1) != 0:
        print(f"Analysis failed at step {step + 1}")
        break

# Storing the coordinates of the original and deflected shape
ox = []  # Original x-coordinates
oy = []  # Original y-coordinates
oz = []  # Original z-coordinates
dx = []  # Deflected x-coordinates
dy = []  # Deflected y-coordinates
dz = []  # Deflected z-coordinates

# Retrieve and store node displacements
displacements = []
for i in range(np + 1):
    disp = ops.nodeDisp(i)
    displacements.append(disp)
    print(f"Node {i} displacement: {disp}")

    # Store original coordinates
    ox.append(ops.nodeCoord(i)[0])
    oy.append(ops.nodeCoord(i)[1])
    oz.append(ops.nodeCoord(i)[2])

    # Store deflected coordinates
    dx.append(ops.nodeCoord(i)[0] + disp[0])
    dy.append(ops.nodeCoord(i)[1] + disp[1])
    dz.append(ops.nodeCoord(i)[2] + disp[2])

# Extract z-direction displacements for 2D plot
z_disp = [disp[2] for disp in displacements]

# Apply exaggeration factor to make deflections more visible
ef = 1e1  # Exaggeration factor

# Plotting the original and deflected shape - 3D Views
fig_3d = plt.figure(figsize=(20, 10))

# 3D Plot with different views
ax_3d = fig_3d.add_subplot(111, projection='3d')
ax_3d.plot(ox, oy, oz, 'o-', label=f'Original {cp.capitalize()} Circle')
ax_3d.plot(dx, dy, [z + ef * dz for z, dz in zip(oz, z_disp)], 'x--', label='Deflected Shape')
ax_3d.set_xlabel('X')
ax_3d.set_ylabel('Y')
ax_3d.set_zlabel('Z')
ax_3d.set_title(f'3D View (Elev=30, Azim=45)')
ax_3d.legend()
ax_3d.grid(True)
ax_3d.view_init(elev=30, azim=45)

plt.suptitle(f'3D Views of Deflection of {cp.capitalize()} Circle')
plt.tight_layout()
plt.show()

# Plotting the original and deflected shape - 2D Views
fig_2d = plt.figure(figsize=(20, 10))

# 2D Plot XY
ax2d_xy = fig_2d.add_subplot(1, 2, 1)
ax2d_xy.plot(ox, oy, 'o-', label='Original Shape')
ax2d_xy.plot(dx, dy, 'x--', label='Deflected Shape')
ax2d_xy.set_xlabel('X')
ax2d_xy.set_ylabel('Y')
ax2d_xy.set_title('2D View in XY Plane')
ax2d_xy.legend()
ax2d_xy.grid(True)

# 2D Plot XZ with exaggerated deflection
ax2d_xz = fig_2d.add_subplot(1, 2, 2)
ax2d_xz.plot(ox, oz, 'o-', label='Original Shape')
ax2d_xz.plot(ox, [z + ef * dz for z, dz in zip(oz, z_disp)], 'x--', label='Deflected Shape (exaggerated)')
ax2d_xz.set_xlabel('X')
ax2d_xz.set_ylabel('Z')
ax2d_xz.set_title('2D View in XZ Plane with Exaggerated Deflection')
ax2d_xz.legend()
ax2d_xz.grid(True)

plt.suptitle(f'2D Views of Deflection of {cp.capitalize()} Circle')
plt.tight_layout()
plt.show()
