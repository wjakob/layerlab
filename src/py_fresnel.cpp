#include <layer/fresnel.h>
#include "python.h"

void python_export_fresnel(py::module &m) {
    /* fresnel.h bindings */
    m.def("fresnelConductor", &fresnelConductor, D(fresnelConductor));
    m.def("fresnelDielectric", [](Float cosThetaI, Float eta) {
        Float cosThetaT;
        Float F = fresnelDielectric(cosThetaI, cosThetaT, eta);
        return std::make_pair(F, cosThetaT);
    }, D(fresnelDielectric));
    m.def("fresnelConductorIntegral", &fresnelConductorIntegral, D(fresnelConductorIntegral));
    m.def("fresnelDielectricIntegral", &fresnelDielectricIntegral, D(fresnelDielectricIntegral));
}
