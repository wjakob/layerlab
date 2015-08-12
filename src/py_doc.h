/*
  This file contains docstrings for the Python bindings.
  Do not edit! These were automatically extracted by mkdoc.py
 */

#define __COUNT(_1, _2, _3, _4, _5, COUNT, ...)  COUNT
#define __VA_SIZE(...)                           __COUNT(__VA_ARGS__, 5, 4, 3, 2, 1)
#define __CAT1(a, b)                             a ## b
#define __CAT2(a, b)                             __CAT1(a, b)
#define __DOC1(n1)                               __doc_##n1
#define __DOC2(n1, n2)                           __doc_##n1##_##n2
#define __DOC3(n1, n2, n3)                       __doc_##n1##_##n2##_##n3
#define __DOC4(n1, n2, n3, n4)                   __doc_##n1##_##n2##_##n3##_##n4
#define __DOC5(n1, n2, n3, n4, n5)               __doc_##n1##_##n2##_##n3##_##n4_##n5
#define DOC(...)                                 __CAT2(__DOC, __VA_SIZE(__VA_ARGS__))(__VA_ARGS__)

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

static const char *__doc_layer_BSDFStorage = R"doc()doc";

static const char *__doc_layer_BSDFStorage_BSDFStorage = R"doc(Map an existing BSDF storage file into memory)doc";

static const char *__doc_layer_BSDFStorage_BSDFStorage_2 = R"doc(Create a new BSDF storage file for the given amount of coefficients etc)doc";

static const char *__doc_layer_BSDFStorage_Header = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_alpha = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_data = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_eta = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_flags = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_identifier = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nBases = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nChannels = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nCoeffs = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nMaxOrder = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nMetadataBytes = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nNodes = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nParameterValues = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_nParameters = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_unused = R"doc()doc";

static const char *__doc_layer_BSDFStorage_Header_version = R"doc()doc";

static const char *__doc_layer_BSDFStorage_alpha = R"doc(Return the Beckmann-equivalent roughness (0: bottom, 1: top surface))doc";

static const char *__doc_layer_BSDFStorage_basisCount = R"doc(Return the number of basis functions stored in this file (usually just 1))doc";

static const char *__doc_layer_BSDFStorage_cdf = 
R"doc(Return a pointer to the coefficients of the CDF associated with the
incident angle ``i``)doc";

static const char *__doc_layer_BSDFStorage_cdf_2 = 
R"doc(Return a pointer to the coefficients of the CDF associated with the
incident angle ``i``)doc";

static const char *__doc_layer_BSDFStorage_channelCount = R"doc(Return the number of color channels)doc";

static const char *__doc_layer_BSDFStorage_close = R"doc(Forcefully release all resources)doc";

static const char *__doc_layer_BSDFStorage_coeff = R"doc(Return the sparse data offset of the given incident and exitant angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffAndCount = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffAndCount_2 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffAndCount_3 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffAndCount_4 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffAndCount_5 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffAndCount_6 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeffCount = R"doc(Return the sparse data size of the given incident and exitant angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeff_2 = R"doc(Return the sparse data offset of the given incident and exitant angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeff_3 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_coeff_4 = 
R"doc(Return the sparse data offset and size of the given incident and exitant
angle pair)doc";

static const char *__doc_layer_BSDFStorage_eta = R"doc(Return the relative index of refraction)doc";

static const char *__doc_layer_BSDFStorage_eval = R"doc(Evaluate the model for the given values of \mu_i, \mu_o, and \phi_d)doc";

static const char *__doc_layer_BSDFStorage_evalLatitudinalAverage = 
R"doc(Evaluate the zeroeth-order Fourier coefficient given size_terpolation
weights)doc";

static const char *__doc_layer_BSDFStorage_evalLatitudinalCDF = 
R"doc(Evaluate the discrete CDF that is used to sample a zenith angle spline
segment)doc";

static const char *__doc_layer_BSDFStorage_extrapolated = 
R"doc(Does this file store coefficients for the harmonic extrapolation-based
model?)doc";

static const char *__doc_layer_BSDFStorage_fromLayer = R"doc(Create a BSDF storage file from a Layer data structure (monochromatic))doc";

static const char *__doc_layer_BSDFStorage_fromLayerGeneral = 
R"doc(Create a BSDF storage file from three Layer data structures (most general
interface))doc";

static const char *__doc_layer_BSDFStorage_fromLayerRGB = R"doc(Create a BSDF storage file from three Layer data structures (RGB))doc";

static const char *__doc_layer_BSDFStorage_getNodes = R"doc(Return the nodes of the underlying discretization in \mu_i and \mu_o)doc";

static const char *__doc_layer_BSDFStorage_interpolateSeries = R"doc(For debugging: return a Fourier series for the given parameters)doc";

static const char *__doc_layer_BSDFStorage_m_cdfMu = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_coeffs = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_filename = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_header = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_metadata = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_mmap = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_nodes = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_offsetTable = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_paramSampleCounts = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_paramSamplePositions = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_paramSamplePositionsNested = R"doc()doc";

static const char *__doc_layer_BSDFStorage_m_reciprocals = R"doc()doc";

static const char *__doc_layer_BSDFStorage_maxOrder = R"doc(Return the number of Fourier coefficients)doc";

static const char *__doc_layer_BSDFStorage_metadata = R"doc(Return metadata attached to the BSDF file (if any))doc";

static const char *__doc_layer_BSDFStorage_nodeCount = R"doc(Return the resolution of the discretization in \mu_i and \mu_o)doc";

static const char *__doc_layer_BSDFStorage_offsetTable = R"doc(Return a posize_ter to the underlying sparse offset table)doc";

static const char *__doc_layer_BSDFStorage_offsetTable_2 = R"doc(Return a posize_ter to the underlying sparse offset table (const version))doc";

static const char *__doc_layer_BSDFStorage_parameterCount = R"doc(Return the number of model parameters)doc";

static const char *__doc_layer_BSDFStorage_parameterSampleCount = R"doc(Return the number of samples associated with parameter ``i``)doc";

static const char *__doc_layer_BSDFStorage_parameterSamplePositions = R"doc(Return the sample positions associated with parameter ``i``)doc";

static const char *__doc_layer_BSDFStorage_pdf = R"doc(Evaluate the model for the given values of \mu_i, \mu_o, and \phi_d)doc";

static const char *__doc_layer_BSDFStorage_sample = R"doc(Importance sample the model)doc";

static const char *__doc_layer_BSDFStorage_setAlpha = R"doc(Return the Beckmann-equivalent roughness (0: bottom, 1: top surface))doc";

static const char *__doc_layer_BSDFStorage_setEta = R"doc(Set the relative index of refraction)doc";

static const char *__doc_layer_BSDFStorage_size = R"doc(Return the size of the underlying representation in bytes)doc";

static const char *__doc_layer_BSDFStorage_stats = R"doc()doc";

static const char *__doc_layer_BSDFStorage_toString = R"doc(Return a string representation)doc";

static const char *__doc_layer_Layer = 
R"doc(Discretized layer reflection model

Describes the linear response to illumination that is incident along a
chosen set of zenith angles. The azimuthal dependence is modeled as an even
real Fourier transform. Each Fourier order is stored in a special
``LayerMode`` data structure.)doc";

static const char *__doc_layer_LayerMode = R"doc(Helper class, which stores one Fourier order of a layer scattering function)doc";

static const char *__doc_layer_LayerMode_LayerMode = 
R"doc(Create a new layer mode data structure

Parameter ``size``:
    Resolution of the underlying matrix, must be a multiple of 2)doc";

static const char *__doc_layer_LayerMode_addScaled = R"doc()doc";

static const char *__doc_layer_LayerMode_clear = R"doc(Reset scattering matrix to that of a clear non-scattering layer)doc";

static const char *__doc_layer_LayerMode_coeff = R"doc(Access a matrix entry)doc";

static const char *__doc_layer_LayerMode_isDiagonal = R"doc(Is this layer represented by a diagonal matrix?)doc";

static const char *__doc_layer_LayerMode_nonZeros = R"doc(Return the number of nonzero coefficients)doc";

static const char *__doc_layer_LayerMode_reflectionBottom = R"doc(Reflection matrix on the bottom interface)doc";

static const char *__doc_layer_LayerMode_reflectionTop = R"doc(Reflection matrix on the top interface)doc";

static const char *__doc_layer_LayerMode_resolution = 
R"doc(Return the underlying matrix resolution (i.e. the number of discretizations
in \mu))doc";

static const char *__doc_layer_LayerMode_reverse = R"doc(Reverse the layer)doc";

static const char *__doc_layer_LayerMode_toString = R"doc(Return a human-readable summary)doc";

static const char *__doc_layer_LayerMode_transmissionBottomTop = R"doc(Transmission matrix going from the bottom to the top interface)doc";

static const char *__doc_layer_LayerMode_transmissionTopBottom = R"doc(Transmission matrix going from the top to the bottom interface)doc";

static const char *__doc_layer_Layer_Layer = 
R"doc(Create a new layer with the given discretization in zenith angles

Parameter ``nodes``:
    A vector with the zenith angle cosines of the chosen discretization

Parameter ``weights``:
    Associated weights for each angle. Usually, the 'nodes' and 'weights'
    are generated using some kind of quadrature rule (e.g. Gauss-Legendre,
    Gauss-Lobatto, etc.)

Parameter ``nFourierOrders``:
    Specifies the number of coefficients to use for the azimuthal Fourier
    expansion (default: 1))doc";

static const char *__doc_layer_Layer_Quartet = R"doc(Helper struct for sparse matrix construction)doc";

static const char *__doc_layer_Layer_Quartet_Quartet = R"doc()doc";

static const char *__doc_layer_Layer_Quartet_i = R"doc()doc";

static const char *__doc_layer_Layer_Quartet_l = R"doc()doc";

static const char *__doc_layer_Layer_Quartet_o = R"doc()doc";

static const char *__doc_layer_Layer_Quartet_value = R"doc()doc";

static const char *__doc_layer_Layer_add = 
R"doc(Combine two layers using the adding equations

Parameter ``l1``:
    Input layer 1

Parameter ``l2``:
    Input layer 2

Parameter ``output``:
    Used to return the resulting layer

Parameter ``homogeneous``:
    When both layers are homogenous, (i.e. if their two sides are
    indistinguishable, this flag should be set to ``True`` to get a speed-
    up). Default: ``False``

Remark:
    In the Python API, the ``output`` parameter is directly returned)doc";

static const char *__doc_layer_Layer_addToBottom = 
R"doc(Append a layer below the current one

This is just a convenience wrapper around Layer::add()

Parameter ``l``:
    The layer to be appended

Parameter ``homogeneous``:
    When the layers are homogenous, (i.e. if their two sides are
    indistinguishable, this flag should be set to ``True`` to get a speed-
    up). Default: ``False``)doc";

static const char *__doc_layer_Layer_addToTop = 
R"doc(Append a layer above the current one

This is just a convenience wrapper around Layer::add()

Parameter ``l``:
    The layer to be appended

Parameter ``homogeneous``:
    When the layers are homogenous, (i.e. if their two sides are
    indistinguishable, this flag should be set to ``True`` to get a speed-
    up). Default: ``False``)doc";

static const char *__doc_layer_Layer_clear = R"doc(Reset scattering matrix to that of a clear non-scattering layer)doc";

static const char *__doc_layer_Layer_eval = R"doc(Evaluate the BSDF for a given pair of zenith angle cosines)doc";

static const char *__doc_layer_Layer_expand = 
R"doc(Solve for the transport matrix of a layer with the given optical thickness
(using Adding-Doubling))doc";

static const char *__doc_layer_Layer_fourierOrders = R"doc(Return the number of Fourier orders)doc";

static const char *__doc_layer_Layer_m_modes = R"doc(Storage for all of the Fourier modes)doc";

static const char *__doc_layer_Layer_m_nodes = R"doc(Integration nodes)doc";

static const char *__doc_layer_Layer_m_weights = R"doc(Integration weights)doc";

static const char *__doc_layer_Layer_matrix = 
R"doc(Return a dense representation of a mode's scattering matrix

All integration weights and cosine foreshortening factors are removed so
that the coefficients can be interpreted as sampled BSDF values)doc";

static const char *__doc_layer_Layer_nodes = R"doc(Return the used integration nodes)doc";

static const char *__doc_layer_Layer_operator_array = R"doc(Look up a mode of the azimuthal Fourier expansion)doc";

static const char *__doc_layer_Layer_operator_array_2 = R"doc(Look up a mode of the azimuthal Fourier expansion (const version))doc";

static const char *__doc_layer_Layer_resolution = R"doc(Return the number of nodes (i.e. the number of discretizations in \mu))doc";

static const char *__doc_layer_Layer_reverse = R"doc(Reverse the layer)doc";

static const char *__doc_layer_Layer_setDiffuse = 
R"doc(Initialize the layer with a diffuse base layer

Parameter ``albedo``:
    The layer's diffuse reflectance (in ``[0, 1]``))doc";

static const char *__doc_layer_Layer_setHenyeyGreenstein = 
R"doc(Initialize the layer with a Henyey-Greenstein phase function

Parameter ``albedo``:
    The layer's single scattering albedo reflectance (in ``[0, 1]``)

Parameter ``g``:
    The layer's HG anisotropy parameter(in ``[-1, 1]``))doc";

static const char *__doc_layer_Layer_setIsotropic = 
R"doc(Initialize the layer with an isotropic phase function

Parameter ``albedo``:
    The layer's single scattering albedo reflectance (in ``[0, 1]``))doc";

static const char *__doc_layer_Layer_setMatusik = 
R"doc(Initialize the layer with a Matusik-style BRDF data file

Parameter ``path``:
    Filename of the BRDF data file

Parameter ``channel``:
    Color channel to extract in ``[0..2]``

Parameter ``fourierOrders``:
    Number of fourier orders that should be used internally in the
    computation. Defaults to the value returned by fourierOrders())doc";

static const char *__doc_layer_Layer_setMicrofacet = 
R"doc(Initialize the layer with a microfacet model (dielectric or conductor)

Parameter ``eta``:
    Relative index of refraction (complex)

Parameter ``alpha``:
    Beckmann roughness coefficient

Parameter ``conserveEnergy``:
    Correct for energy loss due to multiple scattering? Default: ``True``

Parameter ``fourierOrders``:
    Number of fourier orders that should be used internally in the
    computation. Defaults to the value returned by fourierOrders())doc";

static const char *__doc_layer_Layer_setQuartets = R"doc(Initialize from a list of quartets)doc";

static const char *__doc_layer_Layer_setVonMisesFisher = 
R"doc(Initialize the layer with a von Mises-Fisher phase function

Parameter ``albedo``:
    The layer's single scattering albedo reflectance (in ``[0, 1]``)

Parameter ``kappa``:
    The layer's kappa anisotropy parameter)doc";

static const char *__doc_layer_Layer_toString = R"doc(Return a human-readable summary)doc";

static const char *__doc_layer_Layer_weights = R"doc(Return the used integration weights)doc";

static const char *__doc_layer_Layer_write = R"doc(Write the layer coefficients to a sparse file)doc";

static const char *__doc_layer_MemoryMappedFile = R"doc(Basic cross-platform abstraction for memory mapped files)doc";

static const char *__doc_layer_MemoryMappedFile_MemoryMappedFile = R"doc(Create a new memory-mapped file of the specified size)doc";

static const char *__doc_layer_MemoryMappedFile_MemoryMappedFile_2 = R"doc(Map the specified file into memory)doc";

static const char *__doc_layer_MemoryMappedFile_MemoryMappedFile_3 = R"doc(Internal constructor)doc";

static const char *__doc_layer_MemoryMappedFile_createTemporary = 
R"doc(Create a temporary memory-mapped file

Remark:
    When closing the mapping, the file is automatically deleted.)doc";

static const char *__doc_layer_MemoryMappedFile_d = R"doc()doc";

static const char *__doc_layer_MemoryMappedFile_data = R"doc(Return a pointer to the memory-mapped file contents)doc";

static const char *__doc_layer_MemoryMappedFile_data_2 = R"doc(Return a pointer to the memory-mapped file contents (const version))doc";

static const char *__doc_layer_MemoryMappedFile_filename = R"doc(Return the associated filename)doc";

static const char *__doc_layer_MemoryMappedFile_readOnly = R"doc(Return whether the mapped memory region is read-only)doc";

static const char *__doc_layer_MemoryMappedFile_resize = 
R"doc(Resize the memory-mapped file

This involves remapping the file, which will generally change the pointer
obtained via data())doc";

static const char *__doc_layer_MemoryMappedFile_size = R"doc(Return the size of the mapped region)doc";

static const char *__doc_layer_MemoryMappedFile_toString = R"doc(Return a string representation)doc";

static const char *__doc_layer_TColor3 = R"doc(Data type for representing linearized RGB color values)doc";

static const char *__doc_layer_TColor3_TColor3_Scalar_ = R"doc(Initialize the color vector with a uniform value)doc";

static const char *__doc_layer_TColor3_TColor3_Scalar__2 = R"doc(Initialize the color vector with specific per-channel values)doc";

static const char *__doc_layer_TColor3_TColor3_Scalar__3 = R"doc(Assign a color from a dense Eigen expression template)doc";

static const char *__doc_layer_TColor3_b = R"doc(Return a reference to the blue channel)doc";

static const char *__doc_layer_TColor3_b_2 = R"doc(Return a reference to the blue channel (const version))doc";

static const char *__doc_layer_TColor3_clamp = R"doc(Clamp to the positive range)doc";

static const char *__doc_layer_TColor3_g = R"doc(Return a reference to the green channel)doc";

static const char *__doc_layer_TColor3_g_2 = R"doc(Return a reference to the green channel (const version))doc";

static const char *__doc_layer_TColor3_getLuminance = R"doc(Return the associated luminance)doc";

static const char *__doc_layer_TColor3_isValid = R"doc(Check if the color vector contains a NaN/Inf/negative value)doc";

static const char *__doc_layer_TColor3_operator_assign = R"doc(Assign a color from a dense Eigen expression template)doc";

static const char *__doc_layer_TColor3_r = R"doc(Return a reference to the red channel)doc";

static const char *__doc_layer_TColor3_r_2 = R"doc(Return a reference to the red channel (const version))doc";

static const char *__doc_layer_TColor3_toLinearRGB = R"doc(Convert from sRGB to linear RGB)doc";

static const char *__doc_layer_TColor3_toSRGB = R"doc(Convert from linear RGB to sRGB)doc";

static const char *__doc_layer_TFrame3 = 
R"doc(Stores a three-dimensional orthonormal coordinate frame

This class converts between different cartesian coordinate systems and
efficiently computes certain derived quantities based on spherical
coordinates (e.g. cosTheta(), tanTheta(), ..).)doc";

static const char *__doc_layer_TFrame3_TFrame3_Scalar_ = R"doc(Default constructor -- performs no initialization!)doc";

static const char *__doc_layer_TFrame3_TFrame3_Scalar__2 = R"doc(Construct a new coordinate frame from a single vector)doc";

static const char *__doc_layer_TFrame3_TFrame3_Scalar__3 = R"doc(Construct a frame from the given orthonormal vectors)doc";

static const char *__doc_layer_TFrame3_cosPhi = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the cosine of the phi parameter in spherical coordinates)doc";

static const char *__doc_layer_TFrame3_cosPhi2 = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the squared cosine of the phi parameter in spherical coordinates)doc";

static const char *__doc_layer_TFrame3_cosTheta = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the cosine of the angle between the normal and v)doc";

static const char *__doc_layer_TFrame3_cosTheta2 = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the squared cosine of the angle between the normal and v)doc";

static const char *__doc_layer_TFrame3_n = R"doc(Normal direction)doc";

static const char *__doc_layer_TFrame3_operator_eq = R"doc(Equality test)doc";

static const char *__doc_layer_TFrame3_operator_ne = R"doc(Inequality test)doc";

static const char *__doc_layer_TFrame3_s = R"doc(First tangent)doc";

static const char *__doc_layer_TFrame3_sinPhi = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the sine of the phi parameter in spherical coordinates)doc";

static const char *__doc_layer_TFrame3_sinPhi2 = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the squared sine of the phi parameter in spherical coordinates)doc";

static const char *__doc_layer_TFrame3_sinTheta = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the sine of the angle between the normal and v)doc";

static const char *__doc_layer_TFrame3_sinTheta2 = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the squared sine of the angle between the normal and v)doc";

static const char *__doc_layer_TFrame3_t = R"doc(Second tangent)doc";

static const char *__doc_layer_TFrame3_tanTheta = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the tangent of the angle between the normal and v)doc";

static const char *__doc_layer_TFrame3_tanTheta2 = 
R"doc(Assuming that the given direction is in the local coordinate system, return
the squared tangent of the angle between the normal and v)doc";

static const char *__doc_layer_TFrame3_toLocal = R"doc(Convert from world coordinates to local coordinates)doc";

static const char *__doc_layer_TFrame3_toWorld = R"doc(Convert from local coordinates to world coordinates)doc";

static const char *__doc_layer_TNormal3 = R"doc(3-dimensional surface normal representation)doc";

static const char *__doc_layer_TNormal3_TNormal3_Scalar_ = R"doc(Create a new normal with constant component values)doc";

static const char *__doc_layer_TNormal3_TNormal3_Scalar__2 = R"doc(Create a new 3D normal)doc";

static const char *__doc_layer_TNormal3_TNormal3_Scalar__3 = R"doc(Construct a normal from a dense Eigen expression template)doc";

static const char *__doc_layer_TNormal3_operator_assign = R"doc(Assign a normal from a dense Eigen expression template)doc";

static const char *__doc_layer_TPoint = R"doc(Generic N-dimensional point data structure based on Eigen::Matrix)doc";

static const char *__doc_layer_TPoint_TPoint_Scalar__Dimension_ = R"doc(Create a new point with constant component values)doc";

static const char *__doc_layer_TPoint_TPoint_Scalar__Dimension__2 = R"doc(Create a new 2D point (type error if ``Dimension`` != 2))doc";

static const char *__doc_layer_TPoint_TPoint_Scalar__Dimension__3 = R"doc(Create a new 3D point (type error if ``Dimension`` != 3))doc";

static const char *__doc_layer_TPoint_TPoint_Scalar__Dimension__4 = R"doc(Create a new 4D point (type error if ``Dimension`` != 4))doc";

static const char *__doc_layer_TPoint_TPoint_Scalar__Dimension__5 = R"doc(Construct a point from a dense Eigen expression template)doc";

static const char *__doc_layer_TPoint_operator_assign = R"doc(Assign a point from a dense Eigen expression template)doc";

static const char *__doc_layer_TVector = R"doc(Generic N-dimensional vector data structure based on Eigen::Matrix)doc";

static const char *__doc_layer_TVector_TVector_Scalar__Dimension_ = R"doc(Create a new vector with constant component values)doc";

static const char *__doc_layer_TVector_TVector_Scalar__Dimension__2 = R"doc(Create a new 2D vector (type error if ``Dimension`` != 2))doc";

static const char *__doc_layer_TVector_TVector_Scalar__Dimension__3 = R"doc(Create a new 3D vector (type error if ``Dimension`` != 3))doc";

static const char *__doc_layer_TVector_TVector_Scalar__Dimension__4 = R"doc(Create a new 4D vector (type error if ``Dimension`` != 4))doc";

static const char *__doc_layer_TVector_TVector_Scalar__Dimension__5 = R"doc(Construct a vector from a dense Eigen expression template)doc";

static const char *__doc_layer_TVector_operator_assign = R"doc(Assign a vector from a dense Eigen expression template)doc";

static const char *__doc_layer___align_helper = R"doc()doc";

static const char *__doc_layer_convolveFourier = 
R"doc(Computes the Fourier series of a product of even Fourier series using
discrete convolution.

The input series are assumed to have ``ka`` and ``kb`` coefficients, and
the output must have room for ``ka+kb-1`` coefficients.

Remark:
    a First input array of Fourier coefficients

Remark:
    ka Size of the first input array

Remark:
    b Second input array of Fourier coefficients

Remark:
    kb Size of the second input array

Remark:
    c Pointer into the output array

Remark:
    In the Python API, the ``ka`` and ``kb`` parameters are automatically
    computed from the lengths of the input arrays and thus not needed. The
    ``parameter`` is also removed (the result array is directly returned).)doc";

static const char *__doc_layer_coordinateSystem = 
R"doc(Given the unit vector n, find an orthonormal basis {s, t, n}

Based on "Building an Orthonormal Basis from a 3D Unit Vector Without
Normalization" by Jeppe Revall Frisvad in "Journal of Graphics Tools"
16(3), pp. 151-159, August 2012.)doc";

static const char *__doc_layer_evalFourier = 
R"doc(Evaluate an even Fourier series (i.e. containing only cosine terms).

Parameter ``coeffs``:
    Coefficient storage

Parameter ``nCoeffs``:
    Denotes the size of ``coeffs``

Parameter ``phi``:
    Angle for which the series should be evaluated

Returns:
    The value of the series for this angle

Remark:
    In the Python API, the ``nCoeffs`` parameter is automatically computed
    from the length of the input arrays and thus not needed.)doc";

static const char *__doc_layer_evalFourier3 = 
R"doc(Simultaneously evaluate *three* even Fourier series corresponding to the
color channels (Y, R, B) and return a spectral power distribution

Parameter ``coeffs``:
    Coefficient storage

Parameter ``nCoeffs``:
    Denotes the size of ``coeffs``

Parameter ``phi``:
    Angle for which the series should be evaluated

Returns:
    The value of the series for this angle

Remark:
    In the Python API, the ``nCoeffs`` parameter is automatically computed
    from the length of the input arrays and thus not needed.)doc";

static const char *__doc_layer_expCosFourierSeries = 
R"doc(Return Fourier series coefficients for an exponential of a cosine,
specifically the expression "exp(A+B*cos(phi))"

Parameter ``A``:
    The 'A' coefficient in above expression

Parameter ``A``:
    The 'B' coefficient in above expression

Parameter ``relerr``:
    Relative error goal)doc";

static const char *__doc_layer_filonIntegrate = 
R"doc(Compute a Fourier series of the given even function by integrating it
against the basis functions using Filon quadrature

Filon quadrature works by constructing a piecewise quadratic interpolant of
the original function. The Fourier basis functions are then integrated
against this representation, which has an analytic solution. This avoids
all of the problems of traditional quadrature schemes involving highly
oscillatory integrals. It is thus possible to obtain accurate coefficients
even for high orders.

Parameter ``f``:
    Function to be integrated

Parameter ``coeffs``:
    Output buffer used to store the computed coefficients. The function
    adds the computed coefficients to the buffer rather than overwriting
    the existing contents.

Parameter ``nCoeffs``:
    Desired number of coefficients

Parameter ``nEvals``:
    Desired resolution of the piecewise quadratic interpolant

Parameter ``a``:
    Start of the integration, can optionally be set to values other than
    zero. Note that the Fourier basis functions are not orthogonal anymore
    in this case.

Parameter ``b``:
    End of the integration, can be set to values other than pi. Note that
    the Fourier basis functions are not orthogonal anymore in this case.

Remark:
    In the Python API, the ``coeffs`` array is directly returned.)doc";

static const char *__doc_layer_fresnelConductor = 
R"doc(Calculates the unpolarized Fresnel reflection coefficient at a planar
interface having a complex-valued relative index of refraction

Remark:
    The name of this function is a slight misnomer, since it supports the
    general case of a complex-valued relative index of refraction (rather
    than being restricted to conductors)

Parameter ``cosThetaI``:
    Cosine of the angle between the normal and the incident ray

Parameter ``eta``:
    Relative index of refraction (complex))doc";

static const char *__doc_layer_fresnelConductorIntegral = 
R"doc(Calculates the diffuse unpolarized Fresnel reflectance of a conductor

This value quantifies what fraction of diffuse incident illumination will,
on average, be reflected at a conductive material boundary

Parameter ``eta``:
    Relative refractive index (real component)

Parameter ``k``:
    Relative refractive index (imaginary component))doc";

static const char *__doc_layer_fresnelDielectric = 
R"doc(Calculates the unpolarized Fresnel reflection coefficient at a planar
interface between two dielectrics

This function also computes the transmitted direction and returns it using
the ``cosThetaT`` argument. When encountering total internal reflection, it
sets ``cosThetaT=0`` and returns the value 1.

When <tt>cosThetaI < 0</tt>, the function computes the Fresnel reflectance
from the *internal* boundary, which is equivalent to calling the function
with arguments ``fresnelDielectric(abs(cosThetaI), cosThetaT, 1/eta)``.

Remark:
    When accessed from Python, this function has the signature "``F,
    cosThetaT = fresnelDielectric(cosThetaI, eta)``".

Parameter ``cosThetaI``:
    Cosine of the angle between the normal and the incident ray (may be
    negative)

Parameter ``cosThetaT``:
    Argument used to return the cosine of the angle between the normal and
    the transmitted ray, will have the opposite sign of ``cosThetaI``

Parameter ``eta``:
    Relative refractive index)doc";

static const char *__doc_layer_fresnelDielectricIntegral = 
R"doc(Calculates the diffuse unpolarized Fresnel reflectance of a dielectric
material (sometimes referred to as "Fdr").

This value quantifies what fraction of diffuse incident illumination will,
on average, be reflected at a dielectric material boundary

Parameter ``eta``:
    Relative refraction coefficient)doc";

static const char *__doc_layer_fresnelDielectric_2 = 
R"doc(Calculates the unpolarized Fresnel reflection coefficient at a planar
interface between two dielectrics

This is just a convenience wrapper function around the other
``fresnelDielectric`` function, which does not return the transmitted
direction cosine in case it is not needed by the application.

Parameter ``cosThetaI``:
    Cosine of the angle between the normal and the incident ray

Parameter ``eta``:
    Relative refractive index)doc";

static const char *__doc_layer_hg = 
R"doc(Evaluate the HG model using the mu_i, mu_o, phi_d parameterization

Parameter ``mu_i``:
    Incident zenith angle cosine

Parameter ``mu_o``:
    Exitant zenith angle cosine

Parameter ``phi_d``:
    Azimuthal difference angle

Parameter ``g``:
    Anisotropy parameter)doc";

static const char *__doc_layer_hgFourierSeries = 
R"doc(Compute a Fourier series of the HG model

This function first finds the 0-th and 1st-order Fourier coefficients using
elliptic integrals.

The other coefficients can then be determined much more efficiently; the
approach here is based on the idea that the ratio of adjacent coefficients
eventually converges to a constant value. Using a 2nd-order Taylor
expansion, we can obtain a fairly accurate estimate of this ratio somewhere
"in the middle" (i.e. for large $n$, but well before the aforementioned
convergence).

Using a backwards recurrence scheme, we can then determine all previous
ratios and, thereby (using the explicitly computed first values), all
Fourier coefficients.

This approach is based on the article

"A Recurrence Formula For Computing Fourier Components of the Henyey-
Greenstein Phase Function" by E.G. Yanovitskij

Journal of Quantitative Spectroscopy & Radiative Transfer, 57, no 1. 1977

Parameter ``mu_i``:
    Incident zenith angle cosine

Parameter ``mu_o``:
    Exitant zenith angle cosine

Parameter ``g``:
    Anisotropy parameter

Parameter ``kmax``:
    Indicates a desired maximum number of Fourier coefficients. The
    implementation will blur out higher Frequency content to try to achieve
    this number.

Parameter ``result``:
    Storage for the generated Fourier coefficients)doc";

static const char *__doc_layer_math_clamp = R"doc(Generic range clamping function)doc";

static const char *__doc_layer_math_comp_ellint_1 = R"doc(Complete elliptic integral of the first kind (double precision))doc";

static const char *__doc_layer_math_comp_ellint_1_2 = R"doc(Complete elliptic integral of the first kind (single precision))doc";

static const char *__doc_layer_math_comp_ellint_2 = R"doc(Complete elliptic integral of the second kind (double precision))doc";

static const char *__doc_layer_math_comp_ellint_2_2 = R"doc(Complete elliptic integral of the second kind (single precision))doc";

static const char *__doc_layer_math_comp_ellint_3 = R"doc(Complete elliptic integral of the third kind (double precision))doc";

static const char *__doc_layer_math_comp_ellint_3_2 = R"doc(Complete elliptic integral of the third kind (single precision))doc";

static const char *__doc_layer_math_ellint_1 = R"doc(Incomplete elliptic integral of the first kind (double precision))doc";

static const char *__doc_layer_math_ellint_1_2 = R"doc(Incomplete elliptic integral of the first kind (single precision))doc";

static const char *__doc_layer_math_ellint_2 = R"doc(Incomplete elliptic integral of the second kind (double precision))doc";

static const char *__doc_layer_math_ellint_2_2 = R"doc(Incomplete elliptic integral of the second kind (single precision))doc";

static const char *__doc_layer_math_ellint_3 = R"doc(Incomplete elliptic integral of the third kind (double precision))doc";

static const char *__doc_layer_math_ellint_3_2 = R"doc(Incomplete elliptic integral of the first kind (single precision))doc";

static const char *__doc_layer_math_erf = R"doc(Error function (double precision))doc";

static const char *__doc_layer_math_erf_2 = R"doc(Error function (single precision))doc";

static const char *__doc_layer_math_erfinv = R"doc(Inverse error function (double precision))doc";

static const char *__doc_layer_math_erfinv_2 = R"doc(Inverse error function (single precision))doc";

static const char *__doc_layer_math_findInterval = 
R"doc(Find an interval in an ordered set

This function is very similar to ``std::upper_bound``, but it uses a
functor rather than an actual array to permit working with procedurally
defined data. It returns the index ``i`` such that pred(i) is ``True`` and
pred(i+1) is ``False``.

This function is primarily used to locate an interval (i, i+1) for linear
interpolation, hence its name. To avoid issues out of bounds accesses, and
to deal with predicates that evaluate to ``True`` or ``False`` on the
entire domain, the returned left interval index is clamped to the range
``[0, size-2]``.)doc";

static const char *__doc_layer_math_i0e = 
R"doc(Exponentially scaled modified Bessel function of the first kind (order 0),
double precision)doc";

static const char *__doc_layer_math_i0e_2 = 
R"doc(Exponentially scaled modified Bessel function of the first kind (order 0),
single precision)doc";

static const char *__doc_layer_math_legendre_p = R"doc(Evaluate the l-th Legendre polynomial using recurrence, single precision)doc";

static const char *__doc_layer_math_legendre_p_2 = R"doc(Evaluate the l-th Legendre polynomial using recurrence, double precision)doc";

static const char *__doc_layer_math_legendre_p_3 = 
R"doc(Evaluate the an associated Legendre polynomial using recurrence, single
precision)doc";

static const char *__doc_layer_math_legendre_p_4 = 
R"doc(Evaluate the an associated Legendre polynomial using recurrence, double
precision)doc";

static const char *__doc_layer_math_legendre_pd = 
R"doc(Evaluate the l-th Legendre polynomial and its derivative using recurrence,
single precision)doc";

static const char *__doc_layer_math_legendre_pd_2 = 
R"doc(Evaluate the l-th Legendre polynomial and its derivative using recurrence,
double precision)doc";

static const char *__doc_layer_math_legendre_pd_diff = 
R"doc(Evaluate the function ``legendre_pd(l+1, x) - legendre_pd(l-1, x)``, single
precision)doc";

static const char *__doc_layer_math_legendre_pd_diff_2 = 
R"doc(Evaluate the function ``legendre_pd(l+1, x) - legendre_pd(l-1, x)``, double
precision)doc";

static const char *__doc_layer_math_normal_cdf = 
R"doc(Cumulative distribution function of the standard normal distribution
(double precision))doc";

static const char *__doc_layer_math_normal_cdf_2 = 
R"doc(Cumulative distribution function of the standard normal distribution
(single precision))doc";

static const char *__doc_layer_math_normal_quantile = R"doc(Quantile function of the standard normal distribution (double precision))doc";

static const char *__doc_layer_math_normal_quantile_2 = R"doc(Quantile function of the standard normal distribution (single precision))doc";

static const char *__doc_layer_math_safe_acos = 
R"doc(Arccosine variant that gracefully handles arguments > 1 due to roundoff
errors)doc";

static const char *__doc_layer_math_safe_asin = 
R"doc(Arcsine variant that gracefully handles arguments > 1 due to roundoff
errors)doc";

static const char *__doc_layer_math_safe_sqrt = 
R"doc(Square root variant that gracefully handles arguments < 0 due to roundoff
errors)doc";

static const char *__doc_layer_math_signum = 
R"doc(Simple signum function -- note that it returns the FP sign of the input
(and never zero))doc";

static const char *__doc_layer_memString = 
R"doc(Convert a memory amount (in bytes) into a human-readable string
representation

Parameter ``size``:
    An unsigned integer size value

Parameter ``precise``:
    When set to true, a higher-precision string representation is
    generated.)doc";

static const char *__doc_layer_microfacet = 
R"doc(Evaluate the Beckmann distribution-based microfacet BSDF by Walter et al.
using the mu_i, mu_o, phi_d parameterization

Parameter ``mu_i``:
    Incident zenith angle cosine

Parameter ``mu_o``:
    Exitant zenith angle cosine

Parameter ``eta``:
    Relative index of refraction (complex)

Parameter ``alpha``:
    Beckmann roughness parameter

Parameter ``phi_d``:
    Azimuthal difference angle)doc";

static const char *__doc_layer_microfacetFourierSeries = 
R"doc(Compute a Fourier series of the Beckmann-distribution based microfacet BSDF
by Walter et al. (covers both the dielectric and conducting case)

Parameter ``mu_i``:
    Incident zenith angle cosine

Parameter ``mu_o``:
    Exitant zenith angle cosine

Parameter ``eta``:
    Relative index of refraction (complex)

Parameter ``alpha``:
    Beckmann roughness parameter

Parameter ``n``:
    Indicates a desired maximum number of Fourier coefficients. The
    implementation will blur out higher Frequency content to try to achieve
    this number.

Parameter ``relerr``:
    A relative error threshold after which series terms can safely be
    truncated

Parameter ``result``:
    Storage for the generated Fourier coefficients)doc";

static const char *__doc_layer_microfacetNoExp = 
R"doc(Evaluate the Beckmann distribution-based microfacet BSDF by Walter et al.
using the mu_i, mu_o, phi_d parameterization. This version leaves out the
exponential term

Parameter ``mu_i``:
    Incident zenith angle cosine

Parameter ``mu_o``:
    Exitant zenith angle cosine

Parameter ``eta``:
    Relative index of refraction (complex)

Parameter ``k``:
    Absorption coefficient

Parameter ``alpha``:
    Beckmann roughness parameter

Parameter ``phi_d``:
    Azimuthal difference angle)doc";

static const char *__doc_layer_microfacetNoExpFourierSeries = 
R"doc(Compute a Fourier series of the Beckmann-distribution based microfacet BSDF
by Walter et al. (covers both the dielectric and conducting case). This
version leaves out the exponential term

Parameter ``mu_i``:
    Incident zenith angle cosine

Parameter ``mu_o``:
    Exitant zenith angle cosine

Parameter ``eta``:
    Relative index of refraction (complex)

Parameter ``alpha``:
    Beckmann roughness parameter

Parameter ``n``:
    Indicates a desired number of Fourier coefficients (usually 8-12 are
    plenty -- in cases where the function contains high frequencies, it
    will automatically increase ``n``).

Parameter ``phiMax``:
    The implementation minimizes the fitting error on the interval [0,
    phiMax]. If in doubt, set phiMax = math::Pi

Parameter ``result``:
    Storage for the generated Fourier coefficients)doc";

static const char *__doc_layer_pdfFourier = 
R"doc(Evaluate the probability density of the sampling scheme implemented by
sampleFourier() at the position ``phi``.

Parameter ``coeffs``:
    Coefficient storage

Parameter ``nCoeffs``:
    Denotes the size of ``coeffs``

Returns:
    The continuous probability density on [0, \pi])doc";

static const char *__doc_layer_quad_compositeSimpson = 
R"doc(Computes the nodes and weights of a composite Simpson quadrature rule with
the given number of evaluations.

Integration is over the interval $[-1, 1]$, which will be split into $(n-1)
/ 2$ sub-intervals with overlapping endpoints. A 3-point Simpson rule is
applied per interval, which is exact for polynomials of degree three or
less.

Parameter ``n``:
    Desired number of evalution points. Must be an odd number bigger than
    3.

Parameter ``nodes``:
    Length-``n`` array used to store the nodes of the quadrature rule

Parameter ``nodes``:
    Length-``n`` array used to store the weights of the quadrature rule

Remark:
    In the Python API, the ``nodes`` and ``weights`` field are returned as
    a tuple)doc";

static const char *__doc_layer_quad_compositeSimpson38 = 
R"doc(Computes the nodes and weights of a composite Simpson 3/8 quadrature rule
with the given number of evaluations.

Integration is over the interval $[-1, 1]$, which will be split into $(n-1)
/ 3$ sub-intervals with overlapping endpoints. A 4-point Simpson rule is
applied per interval, which is exact for polynomials of degree four or
less.

Parameter ``n``:
    Desired number of evalution points. Must be an odd number bigger than
    3.

Parameter ``nodes``:
    Length-``n`` array used to store the nodes of the quadrature rule

Parameter ``nodes``:
    Length-``n`` array used to store the weights of the quadrature rule

Remark:
    In the Python API, the ``nodes`` and ``weights`` field are returned as
    a tuple)doc";

static const char *__doc_layer_quad_gaussLegendre = 
R"doc(Computes the nodes and weights of a Gauss-Legendre quadrature (aka
"Gaussian quadrature") rule with the given number of evaluations.

Integration is over the interval $[-1, 1]$. Gauss-Legendre quadrature
maximizes the order of exactly integrable polynomials achieves this up to
degree $2n-1$ (where $n$ is the number of function evaluations).

This method is numerically well-behaved until about $n=200$ and then
becomes progressively less accurate. It is generally not a good idea to go
much higher---in any case, a composite or adaptive integration scheme will
be superior for large $n$.

Parameter ``n``:
    Desired number of evalution points

Parameter ``nodes``:
    Length-``n`` array used to store the nodes of the quadrature rule

Parameter ``nodes``:
    Length-``n`` array used to store the weights of the quadrature rule

Remark:
    In the Python API, the ``nodes`` and ``weights`` field are returned as
    a tuple)doc";

static const char *__doc_layer_quad_gaussLobatto = 
R"doc(Computes the nodes and weights of a Gauss-Lobatto quadrature rule with the
given number of evaluations.

Integration is over the interval $[-1, 1]$. Gauss-Lobatto quadrature is
preferable to Gauss-Legendre quadrature whenever the endpoints of the
integration domain should explicitly be included. It maximizes the order of
exactly integrable polynomials subject to this constraint and achieves this
up to degree $2n-3$ (where $n$ is the number of function evaluations).

This method is numerically well-behaved until about $n=200$ and then
becomes progressively less accurate. It is generally not a good idea to go
much higher---in any case, a composite or adaptive integration scheme will
be superior for large $n$.

Parameter ``n``:
    Desired number of evalution points

Parameter ``nodes``:
    Length-``n`` array used to store the nodes of the quadrature rule

Parameter ``nodes``:
    Length-``n`` array used to store the weights of the quadrature rule

Remark:
    In the Python API, the ``nodes`` and ``weights`` field are returned as
    a tuple)doc";

static const char *__doc_layer_sampleFourier = 
R"doc(Sample a angle from an even Fourier series (i.e. containing only cosine
terms).

This is done by importance sampling a uniform piecewise constant
approximation with 2^nRecursions elements, which is constructed on the fly
in a recursive manner.

Parameter ``coeffs``:
    Coefficient storage

Parameter ``nCoeffs``:
    Denotes the size of ``coeffs``

Parameter ``recip``:
    Precomputed array of integer reciprocals, i.e. ``inf, 1, 1/2, 1/3, 1/4,
    ..`` of size ``nCoeffs-1``. This is used to reduce integer-to-FP
    pipeline stalls and division latencies at runtime.

Parameter ``sample``:
    A uniformly distributed random sample in the interval ``[0,1]``

Parameter ``pdf``:
    This parameter is used to return the probability density of the
    sampling scheme evaluated at the generated point on the underlying
    piecewise constant approximation (on [0, \pi])

Parameter ``phi``:
    Used to return the sampled angle (on [0, \pi])

Returns:
    The importance weight (i.e. the value of the Fourier series divided by
    ``pdf``))doc";

static const char *__doc_layer_sampleFourier3 = 
R"doc(Sample a angle from *three* Fourier series corresponding to the color
channels (Y, R, B)

This is done by importance sampling a uniform piecewise constant
approximation with respect to the luminance, which is constructed on the
fly in a recursive manner.

Parameter ``coeffs``:
    Coefficient storage

Parameter ``nCoeffs``:
    Denotes the size of ``coeffs``

Parameter ``recip``:
    Precomputed array of integer reciprocals, i.e. ``inf, 1, 1/2, 1/3, 1/4,
    ..`` of size ``nCoeffs-1``. This is used to reduce integer-to-FP
    pipeline stalls and division latencies at runtime.

Parameter ``sample``:
    A uniformly distributed random sample in the interval ``[0,1]``

Parameter ``pdf``:
    This parameter is used to return the probability density of the
    sampling scheme evaluated at the generated point on the underlying
    piecewise constant approximation (on [0, \pi])

Parameter ``phi``:
    Used to return the sampled angle (on [0, \pi])

Returns:
    The importance weight (i.e. the value of the Fourier series divided by
    ``pdf``))doc";

static const char *__doc_layer_simd_free = R"doc(Free an aligned region of memory)doc";

static const char *__doc_layer_simd_malloc = R"doc(Allocate an aligned region of memory)doc";

static const char *__doc_layer_smithG1 = 
R"doc(Smith's 1D shadowing masking term for the Beckmann microfacet distribution

Parameter ``v``:
    Incident direction

Parameter ``m``:
    Microsurface normal

Parameter ``alpha``:
    Beckmann roughness parameter)doc";

static const char *__doc_layer_spline_eval1D = 
R"doc(Evaluate a cubic spline interpolant of a *uniformly* sampled 1D function

The implementation relies on Catmull-Rom splines, i.e. it uses finite
differences to approximate the derivatives at the endpoints of each spline
segment.

Parameter ``min``:
    Position of the first node

Parameter ``max``:
    Position of the last node

Parameter ``values``:
    Array containing ``size`` regularly spaced evaluations in the range
    [``min``, ``max``] of the approximated function.

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``x``:
    Evaluation point

Parameter ``extrapolate``:
    Extrapolate values when ``x`` is out of range? (default: ``False``)

Remark:
    The Python API lacks the ``size`` parameter, which is inferred
    automatically from the size of the input array.

Remark:
    The Python API provides a vectorized version which evaluates the
    function for many arguments ``x``.

Returns:
    The interpolated value or zero when ``extrapolate=false`` and ``x``
    lies outside of [``min``, ``max``])doc";

static const char *__doc_layer_spline_eval1D_2 = 
R"doc(Evaluate a cubic spline interpolant of a *non*-uniformly sampled 1D
function

The implementation relies on Catmull-Rom splines, i.e. it uses finite
differences to approximate the derivatives at the endpoints of each spline
segment.

Parameter ``nodes``:
    Array containing ``size`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated. They
    must be provided in *increasing* order.

Parameter ``values``:
    Array containing function evaluations matched to the entries of
    ``nodes``.

Parameter ``size``:
    Denotes the size of the ``nodes`` and ``values`` array

Parameter ``x``:
    Evaluation point

Parameter ``extrapolate``:
    Extrapolate values when ``x`` is out of range? (default: ``False``)

Remark:
    The Python API lacks the ``size`` parameter, which is inferred
    automatically from the size of the input array

Remark:
    The Python API provides a vectorized version which evaluates the
    function for many arguments ``x``.

Returns:
    The interpolated value or zero when ``extrapolate=false`` and ``x``
    lies outside of \a [``min``, ``max``])doc";

static const char *__doc_layer_spline_eval2D = 
R"doc(Evaluate a cubic spline interpolant of a uniformly sampled 2D function

This implementation relies on a tensor product of Catmull-Rom splines, i.e.
it uses finite differences to approximate the derivatives for each
dimension at the endpoints of spline patches.

Parameter ``nodes1``:
    Arrays containing ``size1`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated on
    the ``X`` axis (in increasing order)

Parameter ``size1``:
    Denotes the size of the ``nodes1`` array

Parameter ``nodes``:
    Arrays containing ``size2`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated on
    the ``Y`` axis (in increasing order)

Parameter ``size2``:
    Denotes the size of the ``nodes2`` array

Parameter ``values``:
    A 2D floating point array of ``size1*size2`` cells containing
    irregularly spaced evaluations of the function to be interpolated.
    Consecutive entries of this array correspond to increments in the ``X``
    coordinate.

Parameter ``x``:
    ``X`` coordinate of the evaluation point

Parameter ``y``:
    ``Y`` coordinate of the evaluation point

Parameter ``extrapolate``:
    Extrapolate values when ``p`` is out of range? (default: ``False``)

Remark:
    The Python API lacks the ``size1`` and ``size2`` parameters, which are
    inferred automatically from the size of the input arrays.

Returns:
    The interpolated value or zero when ``extrapolate=false``tt> and
    ``(x,y)`` lies outside of the node range)doc";

static const char *__doc_layer_spline_evalSpline = 
R"doc(Compute the definite integral and derivative of a cubic spline that is
parameterized by the function values and derivatives at the endpoints of
the interval ``[0, 1]``.

Parameter ``f0``:
    The function value at the left position

Parameter ``f1``:
    The function value at the right position

Parameter ``d0``:
    The function derivative at the left position

Parameter ``d1``:
    The function derivative at the right position

Parameter ``t``:
    The parameter variable

Returns:
    The interpolated function value at ``t``)doc";

static const char *__doc_layer_spline_evalSplineD = 
R"doc(Compute the value and derivative of a cubic spline that is parameterized by
the function values and derivatives of the interval ``[0, 1]``.

Parameter ``f0``:
    The function value at the left position

Parameter ``f1``:
    The function value at the right position

Parameter ``d0``:
    The function derivative at the left position

Parameter ``d1``:
    The function derivative at the right position

Parameter ``t``:
    The parameter variable

Returns:
    The interpolated function value and its derivative at ``t``)doc";

static const char *__doc_layer_spline_evalSplineI = 
R"doc(Compute the definite integral and value of a cubic spline that is
parameterized by the function values and derivatives of the interval ``[0,
1]``.

Parameter ``f0``:
    The function value at the left position

Parameter ``f1``:
    The function value at the right position

Parameter ``d0``:
    The function derivative at the left position

Parameter ``d1``:
    The function derivative at the right position

Returns:
    The definite integral and the interpolated function value at ``t``)doc";

static const char *__doc_layer_spline_evalSplineWeights = 
R"doc(Compute weights to perform a spline-interpolated lookup on a *uniformly*
sampled 1D function.

The implementation relies on Catmull-Rom splines, i.e. it uses finite
differences to approximate the derivatives at the endpoints of each spline
segment. The resulting weights are identical those internally used by
sample1D().

Parameter ``min``:
    Position of the first node

Parameter ``max``:
    Position of the last node

Parameter ``size``:
    Denotes the number of function samples

Parameter ``x``:
    Evaluation point

Parameter ``weights``:
    Pointer to a weight array of size 4 that will be populated

Parameter ``offset``:
    Offset into the function samples associated with weights[0]

Parameter ``extrapolate``:
    Extrapolate values when ``x`` is out of range? (default: ``False``)

Remark:
    In the Python API, the ``offset`` and ``weights`` parameters are
    returned as the second and third elements of a triple.

Returns:
    ``True`` on success and ``False`` when ``extrapolate=false`` and ``x``
    lies outside of [``min``, ``max``])doc";

static const char *__doc_layer_spline_evalSplineWeights_2 = 
R"doc(Compute weights to perform a spline-interpolated lookup on a
*non*-uniformly sampled 1D function.

The implementation relies on Catmull-Rom splines, i.e. it uses finite
differences to approximate the derivatives at the endpoints of each spline
segment. The resulting weights are identical those internally used by
sample1D().

Parameter ``nodes``:
    Array containing ``size`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated. They
    must be provided in *increasing* order.

Parameter ``size``:
    Denotes the size of the ``nodes`` array

Parameter ``x``:
    Evaluation point

Parameter ``weights``:
    Pointer to a weight array of size 4 that will be populated

Parameter ``offset``:
    Offset into the function samples associated with weights[0]

Parameter ``extrapolate``:
    Extrapolate values when ``x`` is out of range? (default: ``False``)

Remark:
    The Python API lacks the ``size`` parameter, which is inferred
    automatically from the size of the input array. The ``offset`` and
    ``weights`` parameters are returned as the second and third elements of
    a triple.

Returns:
    ``True`` on success and ``False`` when ``extrapolate=false`` and ``x``
    lies outside of [``min``, ``max``])doc";

static const char *__doc_layer_spline_integrate1D = 
R"doc(Computes a prefix sum of integrals over segments of a *uniformly* sampled
1D Catmull-Rom spline interpolant

This is useful for sampling spline segments as part of an importance
sampling scheme (in conjunction with sample1D)

Parameter ``min``:
    Position of the first node

Parameter ``max``:
    Position of the last node

Parameter ``values``:
    Array containing ``size`` regularly spaced evaluations in the range
    [``min``, ``max``] of the approximated function.

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``out``:
    An array with ``size`` entries, which will be used to store the prefix
    sum

Remark:
    The Python API lacks the ``size`` and ``out`` parameters. The former is
    inferred automatically from the size of the input array, and ``out`` is
    returned as a list.)doc";

static const char *__doc_layer_spline_integrate1D_2 = 
R"doc(Computes a prefix sum of integrals over segments of a *non*-uniformly
sampled 1D Catmull-Rom spline interpolant

This is useful for sampling spline segments as part of an importance
sampling scheme (in conjunction with sample1D)

Parameter ``nodes``:
    Array containing ``size`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated. They
    must be provided in *increasing* order.

Parameter ``values``:
    Array containing function evaluations matched to the entries of
    ``nodes``.

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``out``:
    An array with ``size`` entries, which will be used to store the prefix
    sum

Remark:
    The Python API lacks the ``size`` and ``out`` parameters. The former is
    inferred automatically from the size of the input array, and ``out`` is
    returned as a list.)doc";

static const char *__doc_layer_spline_invert1D = 
R"doc(Invert a cubic spline interpolant of a *uniformly* sampled 1D function. The
spline interpolant must be *monotonically increasing*.

Parameter ``min``:
    Position of the first node

Parameter ``max``:
    Position of the last node

Parameter ``values``:
    Array containing ``size`` regularly spaced evaluations in the range
    [``min``, ``max``] of the approximated function.

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``y``:
    Input parameter for the inversion

Returns:
    The spline parameter ``t`` such that ``eval1D(..., t)=y``)doc";

static const char *__doc_layer_spline_invert1D_2 = 
R"doc(Invert a cubic spline interpolant of a *non*-uniformly sampled 1D function.
The spline interpolant must be *monotonically increasing*.

Parameter ``nodes``:
    Array containing ``size`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated. They
    must be provided in *increasing* order.

Parameter ``values``:
    Array containing function evaluations matched to the entries of
    ``nodes``.

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``y``:
    Input parameter for the inversion

Returns:
    The spline parameter ``t`` such that ``eval1D(..., t)=y``)doc";

static const char *__doc_layer_spline_sample1D = 
R"doc(Importance sample a segment of a *uniformly* sampled 1D Catmull-Rom spline
interpolant

Parameter ``min``:
    Position of the first node

Parameter ``max``:
    Position of the last node

Parameter ``values``:
    Array containing ``size`` regularly spaced evaluations in the range
    [``min``, ``max``] of the approximated function.

Parameter ``cdf``:
    Array containing a cumulative distribution function computed by
    integrate1D().

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``sample``:
    A uniformly distributed random sample in the interval ``[0,1]``

Parameter ``fval``:
    If set to a non-null pointer, this argument will be used to return the
    value of the spline at the sampled position

Parameter ``pdf``:
    If set to a non-null pointer, this argument will be used to return the
    probability density at the sampled position (which only differs from
    ``fval`` when the function does not integrate to 1)

Remark:
    The Python API lacks the ``size``, ``fval`` and ``pdf`` parameters. The
    first is automatically inferred from the size of the input array, and
    ``fval`` and ``pdf`` are returned as the second and third element of
    the return value, which is now a tuple.

Returns:
    The sampled position)doc";

static const char *__doc_layer_spline_sample1D_2 = 
R"doc(Importance sample a segment of a *non*-uniformly sampled 1D Catmull-Rom
spline interpolant

Parameter ``nodes``:
    Array containing ``size`` non-uniformly spaced values denoting
    positions the where the function to be interpolated was evaluated. They
    must be provided in *increasing* order.

Parameter ``values``:
    Array containing function evaluations matched to the entries of
    ``nodes``.

Parameter ``cdf``:
    Array containing a cumulative distribution function computed by
    integrate1D().

Parameter ``size``:
    Denotes the size of the ``values`` array

Parameter ``sample``:
    A uniformly distributed random sample in the interval ``[0,1]``

Parameter ``fval``:
    If set to a non-null pointer, this argument will be used to return the
    value of the spline at the sampled position

Parameter ``pdf``:
    If set to a non-null pointer, this argument will be used to return the
    probability density at the sampled position (which only differs from
    ``fval`` when the function does not integrate to 1)

Remark:
    The Python API lacks the ``size``, ``fval`` and ``pdf`` parameters. The
    first is automatically inferred from the size of the input array, and
    ``fval`` and ``pdf`` are returned as the second and third element of
    the return value, which is now a tuple.

Returns:
    The sampled position)doc";

static const char *__doc_layer_timeString = 
R"doc(Convert a time difference (in seconds) into a human-readable string
representation

Parameter ``time``:
    Time difference in (fractional) sections

Parameter ``precise``:
    When set to true, a higher-precision string representation is
    generated.)doc";
