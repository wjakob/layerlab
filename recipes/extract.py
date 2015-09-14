import jsonpickle
name = "Cr"
l1 = open("%s.eta.spd" % name).readlines()
l2 = open("%s.k.spd" % name).readlines()

values = []
wavelengths = []

for i in range(len(l1)):
    if l1[i].startswith('#'):
        continue
    lambda1, x1 = l1[i].strip().split()
    lambda2, x2 = l2[i].strip().split()
    if lambda1 != lambda2:
        raise Exception("internal error")
    values.append(complex(x1) + complex(x2)*1j)
    wavelengths.append(float(lambda1))
print(wavelengths)
print(values)
