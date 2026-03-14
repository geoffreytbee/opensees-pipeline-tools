# Pipe-Soil Interaction Theory

## What is Pipe-Soil Interaction?

Pipe-soil interaction (PSI) is the mechanical coupling between a buried pipeline and the surrounding soil. When a pipeline moves relative to the soil — whether due to ground displacement, thermal expansion, or internal pressure — the soil resists that movement. The magnitude and distribution of soil resistance govern the pipeline's structural response: its displacements, internal forces, and stresses.

Accurate modeling of PSI is essential for buried pipeline design and integrity assessment because:

1. **The soil is the load path.** Unlike above-ground structures, buried pipelines receive their loads through the soil. Ground movement from faults, landslides, or liquefaction is transmitted to the pipe through soil resistance.
2. **The response is nonlinear.** Soil resistance is displacement-dependent and reaches an ultimate value beyond which the soil yields. This nonlinearity must be captured to predict peak pipe stresses accurately.
3. **The interaction is directional.** Soil resistance differs in each direction (axial, lateral, vertical bearing, uplift) due to different failure mechanisms in the soil.

## The Four Spring Directions

In the discrete spring approach used by this toolkit, the pipe-soil interface is represented by independent nonlinear springs in four directions at each node along the pipeline:

### Axial Springs (t-z)

Axial springs resist relative movement between the pipe and soil along the pipe axis. The resistance comes from friction on the pipe-soil interface and depends on:

- Pipe diameter and coating (which affects the interface friction angle)
- Burial depth (which determines the normal stress on the pipe surface)
- Soil friction angle and cohesion

The ultimate axial force per unit length is given by ALA 2001 Section 7.2:

$$t_u = \pi D \left( \frac{1}{2} \gamma H (1 + K_0) \tan \delta + c_a \right)$$

where $D$ is the pipe outside diameter, $\gamma$ is the soil unit weight, $H$ is the depth to the pipe centerline, $K_0 = 1 - \sin\phi$ is the at-rest earth pressure coefficient, $\delta = f \cdot \phi$ is the interface friction angle, and $c_a$ is the pipe-soil adhesion.

### Lateral Horizontal Springs (p-y)

Lateral springs resist pipe movement perpendicular to its axis in the horizontal plane. The failure mechanism is analogous to a laterally loaded pile, with passive earth pressure developing in front of the pipe and active (or no) pressure behind it.

The ultimate lateral force per unit length (ALA 2001 Section 7.3):

$$p_u = N_{ch} \cdot c \cdot D + N_{qh} \cdot \gamma \cdot H \cdot D$$

where $N_{ch}$ and $N_{qh}$ are bearing capacity factors that depend on the embedment ratio $H/D$ and friction angle $\phi$.

### Vertical Bearing Springs (q-z downward)

Bearing springs resist downward pipe movement into the soil. The ultimate capacity follows the general bearing capacity equation adapted for a strip footing at depth (ALA 2001 Section 7.4):

$$q_u = N_c \cdot c \cdot D + N_q \cdot \gamma \cdot H \cdot D + \frac{1}{2} N_\gamma \cdot \gamma \cdot D^2$$

where $N_c$, $N_q$, and $N_\gamma$ are Terzaghi/Meyerhof bearing capacity factors.

### Uplift Springs (vertical upward)

Uplift springs resist upward pipe movement out of the soil. The resistance comes from the weight of the soil wedge above the pipe and any cohesive resistance along the failure surface (ALA 2001 Section 7.5):

$$Q_u = N_{cv} \cdot c \cdot D + N_{qv} \cdot \gamma \cdot H \cdot D$$

The uplift capacity is capped at the total weight of soil above the pipe ($\gamma H D$) since the failure mechanism cannot mobilize more resistance than the overburden.

## The ALA 2001 Bilinear Backbone Curve

The ALA 2001 guidelines model each spring as **bilinear elastic-perfectly-plastic**:

1. **Elastic region:** The spring force increases linearly with displacement up to the yield displacement. The initial stiffness is $k = F_u / \delta_y$, where $F_u$ is the ultimate force and $\delta_y$ is the yield displacement.

2. **Plastic region:** Beyond the yield displacement, the spring force remains constant at the ultimate value regardless of further displacement. The soil has yielded and can provide no additional resistance.

```
Force
  ^
  |         _______________
  |        /
  |       /
  |      /  <- elastic region
  |     /
  |    /
  |   /
  +---+---+---+---+---+---> Displacement
      δy
```

This simple model captures the essential nonlinear behavior of soil resistance while requiring only two parameters per spring (ultimate force and yield displacement). More sophisticated multi-linear or hyperbolic models exist but the bilinear model is standard practice and is the foundation of this toolkit.

## Tributary Length Scaling

A critical step in applying soil spring properties to a finite element model is **tributary length scaling**. The ALA 2001 equations give spring properties per unit length of pipe (e.g., $t_u$ in kN/m). However, in the finite element model, springs are attached at discrete nodes, not continuously along the pipe.

Each node's spring must represent the resistance over its **tributary length** — the length of pipe "served" by that node:

- **Interior nodes:** tributary length = average of the two adjacent element lengths
- **End nodes:** tributary length = half the adjacent element length

The per-node spring force is then:

$$F_{u,\text{node}} = t_u \times L_{\text{trib}}$$

And the per-node spring stiffness:

$$k_{\text{node}} = \frac{F_{u,\text{node}}}{\delta_y} = \frac{t_u \times L_{\text{trib}}}{\delta_y}$$

!!! warning "Critical Convention"
    All spring stiffness values assigned to OpenSees zero-length elements must be **per node** (already scaled by tributary length). The `springs.py` module returns per-unit-length values. Scaling happens during model assembly. Getting this wrong silently produces incorrect results.

The sum of all tributary lengths must equal the total pipe length. This ensures that the total soil resistance in the discrete model equals the integral of the continuous spring resistance over the pipe length — a necessary condition for the discrete model to converge to the continuous solution as the mesh is refined.

## Abaqus PSI Elements vs. OpenSees Zero-Length Springs

Commercial FEA codes like Abaqus provide dedicated pipe-soil interaction (PSI) elements that automatically handle the spring distribution and tributary length scaling. In Abaqus, PSI elements are defined between a pipe node and a ground node, and the element formulation internally accounts for the element length in the spring stiffness calculation.

In OpenSees, there is no dedicated PSI element. Instead, zero-length elements with appropriate material models (e.g., `ElasticPP` for bilinear springs) are used. The user is responsible for:

1. Creating a fixed "ground" node at the same location as each pipe node
2. Attaching a zero-length element between the pipe node and ground node
3. Assigning the spring material with stiffness already scaled by tributary length
4. Constraining the ground node in all degrees of freedom

For small strains and displacements, the two approaches produce equivalent results. The Benchmark 01 validation demonstrates agreement within 3.1% for bending stress, which is within the expected range for different element formulations and numerical implementations.

The OpenSees approach offers more flexibility — any constitutive model can be used for the spring material, and the springs can be modified independently at each node (e.g., for spatially varying soil properties or probabilistic analysis). This flexibility is a key advantage for research and advanced engineering applications.
