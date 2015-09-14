# Creates a rough gold layer with a rough dielectric coating containing an
# anisotropic scattering medium

import sys
sys.path.append('.')

from utils.materials import gold
from utils.cie import get_rgb

import layerlab as ll

eta_top = 1.5

# This step integrates the spectral IOR against the CIE XYZ curves to obtain
# equivalent sRGB values. This may seem fairly approximate but turns out to
# yield excellent agreement with spectral reference renders
print('Computing gold IOR parameters')
eta_bot = get_rgb(gold)

alpha_top = 0.1  # Beckmann roughness of top layer (coating)
alpha_bot = 0.1  # Beckmann roughness of bottom layer (gold)

# Medium parameters
g = 0.5  # Scattering anisotropy
albedo = [0.25, 0.0, 0.95]  # Single scattering albedo
tau = 0.5  # Optical depth

# Construct quadrature scheme suitable for the material
n_top, m_top = ll.parameterHeuristicMicrofacet(eta=eta_top, alpha=alpha_top)
n_bot, m_bot = ll.parameterHeuristicMicrofacet(eta=eta_bot[0], alpha=alpha_bot)
n_med, m_med = ll.parameterHeuristicHG(g=g)

n = max(n_top, n_bot)  # Max of zenith angle discretization
m = m_top              # Number of Fourier orders determined by top layer
mu, w = ll.quad.gaussLobatto(n)
print("# of nodes = %i, fourier orders = %i" % (n, m))

# Construct coating layer
print("Creating coating layer")
coating = ll.Layer(mu, w, m)
coating.setMicrofacet(eta=eta_top, alpha=alpha_top)

output = []
for channel in range(3):
    # Construct diffuse bottom layer for each channel
    print("Creating metal layer")
    l = ll.Layer(mu, w, m)
    l.setMicrofacet(eta=eta_bot[channel], alpha=alpha_bot)

    # Construct medium layer
    print("Creating medium layer")
    l2 = ll.Layer(mu, w, m)
    l2.setHenyeyGreenstein(g=g, albedo=albedo[channel])
    l2.expand(tau)

    # Apply medium layer
    print("Applying medium ..")
    l.addToTop(l2)

    # Apply coating
    print("Applying coating..")
    l.addToTop(coating)
    output.append(l)

# .. and write to disk
print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("output.bsdf", *output)
storage.close()
