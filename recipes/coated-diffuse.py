import sys
sys.path.append('.')

import layerlab as ll

albedo = [0.5, 0.5, 0.6]
eta    = 1.5
alpha  = 0.1 # Beckmann roughness

# Construct quadrature scheme suitable for the material
n, m = ll.parameterHeuristicMicrofacet(eta = eta, alpha = alpha)
mu, w = ll.quad.gaussLobatto(n)
print("# of nodes = %i, fourier orders = %i" % (n, m))

# Construct coating layer
print("Creating coating layer")
coating = ll.Layer(mu, w, m)
coating.setMicrofacet(eta = eta, alpha = alpha)

output = []
for channel in range(3):
    # Construct diffuse bottom layer for each channel
    print("Creating diffuse layer")
    l = ll.Layer(mu, w, m)
    l.setDiffuse(albedo[channel])
    
    # Apply coating
    print("Applying coating..")
    l.addToTop(coating)
    output.append(l)

# .. and write to disk
print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("output.bsdf", output[0], output[1], output[2])
storage.close()
