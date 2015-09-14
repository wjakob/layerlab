#include <layer/fresnel.h>
#include "python.h"

void python_export_fresnel(py::module &m) {
    /* fresnel.h bindings */
    m.def("fresnelConductor", &fresnelConductor, D(fresnelConductor),
          py::arg("cosThetaI"), py::arg("eta"));
    m.def("fresnelDielectric", [](Float cosThetaI, Float eta) {
        Float cosThetaT;
        Float F = fresnelDielectric(cosThetaI, cosThetaT, eta);
        return std::make_pair(F, cosThetaT);
    }, D(fresnelDielectric), py::arg("cosThetaI"), py::arg("eta"));
    m.def("fresnelConductorIntegral", &fresnelConductorIntegral, D(fresnelConductorIntegral));
    m.def("fresnelDielectricIntegral", &fresnelDielectricIntegral, D(fresnelDielectricIntegral));
}
