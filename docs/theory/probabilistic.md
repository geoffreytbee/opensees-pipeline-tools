# Probabilistic Methods

!!! info "Planned Feature"
    The Monte Carlo simulation framework described on this page is planned but not yet implemented. This page documents the theoretical foundation and intended approach.

## Why Probabilistic Analysis?

Deterministic pipeline analysis uses single "best estimate" or "design" values for soil properties, pipe wall thickness, and ground displacement profiles. In reality, these quantities are uncertain and spatially variable:

- **Soil properties** (friction angle, unit weight, cohesion) vary across a site and are known only at discrete borehole locations
- **Pipe wall thickness** varies due to manufacturing tolerances and corrosion
- **Ground displacement profiles** from fault rupture or landslide movement are inherently uncertain

Probabilistic analysis using Monte Carlo simulation accounts for these uncertainties by running the analysis many times with different input values sampled from their probability distributions. The result is not a single stress value but a distribution of stress values, from which engineers can assess the probability of exceeding a given limit state.

## Random Fields

Many soil properties are not just uncertain — they are **spatially correlated**. The friction angle at one point along the pipeline is likely to be similar to the friction angle at a nearby point, but may be quite different from the friction angle at a distant point. This spatial structure is captured by **random fields**.

A random field is characterized by:

1. **Marginal distribution** — the probability distribution at any single point (e.g., normal with mean 35 degrees and standard deviation 3 degrees for friction angle)
2. **Correlation function** — how the property values at two points are related as a function of the distance between them

The key parameter of the correlation function is the **spatial correlation length** (also called the **scale of fluctuation**). This is the distance over which property values remain significantly correlated:

- **Short correlation length** (e.g., 5 m): properties vary rapidly along the pipeline, averaging out over the pipe length
- **Long correlation length** (e.g., 500 m): properties vary slowly, meaning long sections of pipe may all experience above- or below-average conditions simultaneously

The correlation length has a profound effect on pipeline response. Very short correlation lengths lead to spatial averaging that reduces the effect of variability. Very long correlation lengths approach the deterministic case with a single uncertain value applied everywhere. Intermediate correlation lengths — typically on the order of the pipeline feature length (fault crossing width, landslide extent) — often produce the most critical response.

## Planned Implementation

The planned Monte Carlo framework will support two methods for generating spatially correlated random fields:

### Cholesky Decomposition

The most straightforward approach:

1. Define the correlation matrix $\mathbf{C}$ where $C_{ij} = \rho(|x_i - x_j|)$ and $\rho$ is the correlation function
2. Compute the Cholesky decomposition $\mathbf{C} = \mathbf{L}\mathbf{L}^T$
3. Generate a vector of independent standard normal random variables $\mathbf{z}$
4. Compute the correlated field as $\mathbf{y} = \boldsymbol{\mu} + \mathbf{L}\mathbf{z} \cdot \boldsymbol{\sigma}$

This method is exact but has $O(n^3)$ computational cost for $n$ nodes, which may be prohibitive for very large meshes.

### Karhunen-Loeve Expansion

A spectral method that represents the random field as a series:

$$H(x, \theta) = \mu(x) + \sum_{i=1}^{M} \sqrt{\lambda_i} \, \phi_i(x) \, \xi_i(\theta)$$

where $\lambda_i$ and $\phi_i$ are eigenvalues and eigenfunctions of the correlation function, and $\xi_i$ are uncorrelated standard normal random variables. The series is truncated at $M$ terms, providing a reduced-order representation that is computationally efficient.

### Target Applications

- Soil spring property variability (friction angle, unit weight, cohesion)
- Pipe wall thickness variability along the pipeline
- Ground displacement profile uncertainty (amplitude, width, location)
- Combined uncertainty in multiple parameters

The framework will collect results from all Monte Carlo realizations and provide statistical summaries (mean, standard deviation, percentiles) and cumulative distribution functions for key response quantities such as peak bending stress and maximum displacement.
