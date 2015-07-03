#include <layer/vector.h>
#include <layer/frame.h>
#include <layer/color.h>
#include <Eigen/Dense>

#include "python.h"

template <typename Type> void init_fixed_from_buffer(Type &v, py::buffer &b) {
    typedef typename Type::Scalar Scalar;

    py::buffer_info info = b.request();
    if (info.format != py::format_descriptor<Scalar>::value())
        throw std::runtime_error("Incompatible buffer format!");
    if (!((info.ndim == 1 && info.strides[0] == sizeof(Scalar)) ||
          (info.ndim == 2 &&
              ((info.shape[0] == 1 && info.strides[0] == sizeof(Scalar) &&
                info.shape[1] == Type::Dimension) ||
               (info.shape[1] == 1 && info.strides[1] == sizeof(Scalar) &&
                info.shape[0] == Type::Dimension)))))
        throw std::runtime_error("Incompatible buffer dimension!");

    memcpy(v.data(), info.ptr, sizeof(Scalar) * Type::Dimension);
}

/// Creates Python bindings for an Eigen order-1 tensor of size 3 (i.e. a vector/normal/point)
template <typename Type>
py::class_<Type> bind_eigen_1_3(py::module &m, const char *name,
                                py::object parent = py::object()) {
    typedef typename Type::Scalar Scalar;

    py::class_<Type> vector(m, name, parent);
    vector
        /* Constructors */
        .def(py::init<>())
        .def(py::init<Scalar>())
        .def(py::init<Scalar, Scalar, Scalar>())
        .def("__init__", [](Type &v, const std::vector<Scalar> &v2) {
            if (v2.size() != Type::Dimension)
                throw std::runtime_error("Incompatible size!");
            memcpy(v.data(), &v2[0], sizeof(Scalar) * Type::Dimension);
        })
        .def("__init__", [](Type &v, py::buffer b) {
            init_fixed_from_buffer(v, b);
        })

        /* Initialization */
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })
        .def("setZero", [](Type &m) { m.setZero(); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Matrix to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * Scalar())
        .def_cast(py::self / Scalar())
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)

        /* Python protocol implementations */
        .def("__len__", [](const Type &) { return (int) Type::Dimension; })
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v;
            return oss.str();
        })
        .def("__getitem__", [](const Type &c, int i) {
            if (i < 0 || i >= Type::Dimension)
                throw py::index_error();
            return c[i];
         })
        .def("__setitem__", [](Type &c, int i, Scalar v) {
             if (i < 0 || i >= Type::Dimension)
                 throw py::index_error();
            c[i] = v;
         })

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),        /* Pointer to buffer */
                sizeof(Scalar),  /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::value(),
                1, { (size_t) Type::Dimension },
                { sizeof(Scalar) }
            );
        });
    return vector;
}

/// Creates Python bindings for a dynamic Eigen order-1 tensor (i.e. a vector)
template <typename Type>
py::class_<Type> bind_eigen_1(py::module &m, const char *name,
                              py::object parent = py::object()) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> vector(m, name, parent);

    vector
        /* Constructors */
        .def(py::init<>())
        .def(py::init<size_t>())
        .def("__init__", [](Type &v, const std::vector<Scalar> &v2) {
            new (&v) Type(v2.size());
            memcpy(v.data(), &v2[0], sizeof(Scalar) * v2.size());
        })
        .def("__init__", [](Type &v, py::buffer b) {
            py::buffer_info info = b.request();
            if (info.format != py::format_descriptor<Scalar>::value()) {
                throw std::runtime_error("Incompatible buffer format!");
            } else if (info.ndim == 1 && info.strides[0] == sizeof(Scalar)) {
                new (&v) Type(info.shape[0]);
                memcpy(v.data(), info.ptr, sizeof(Scalar) * info.shape[0]);
            } else if (info.ndim == 2 && ((info.shape[0] == 1 && info.strides[0] == sizeof(Scalar))
                                       || (info.shape[1] == 1 && info.strides[1] == sizeof(Scalar)))) {
                new (&v) Type(info.shape[0] * info.shape[1]);
                memcpy(v.data(), info.ptr, sizeof(Scalar) * info.shape[0] * info.shape[1]);
            } else {
                throw std::runtime_error("Incompatible buffer dimension!");
            }
        })

        /* Size query functions */
        .def("size", &Type::size)
        .def("cols", &Type::cols)
        .def("rows", &Type::rows)

        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })

        /* Resizing */
        .def("resize", [](Type &m, size_t s0) { m.resize(s0); })
        .def("resizeLike", [](Type &m, const Type &m2) { m.resizeLike(m2); })
        .def("conservativeResize", [](Type &m, size_t s0) { m.conservativeResize(s0); })

        /* Component-wise operations */
        .def("cwiseAbs", &Type::cwiseAbs)
        .def("cwiseAbs2", &Type::cwiseAbs2)
        .def("cwiseSqrt", &Type::cwiseSqrt)
        .def("cwiseInverse", &Type::cwiseInverse)
        .def("cwiseMin", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMin(m2); })
        .def("cwiseMax", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMax(m2); })
        .def("cwiseMin", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMin(s); })
        .def("cwiseMax", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMax(s); })
        .def("cwiseProduct", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseProduct(m2); })
        .def("cwiseQuotient", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseQuotient(m2); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * Scalar())
        .def_cast(py::self / Scalar())

        /* Arithmetic in-place operators */
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)

        /* Python protocol implementations */
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v.transpose();
            return oss.str();
        })
        .def("__getitem__", [](const Type &m, size_t i) {
            if (i >= (size_t) m.size())
                throw py::index_error();
            return m[i];
         })
        .def("__setitem__", [](Type &m, size_t i, Scalar v) {
            if (i >= (size_t) m.size())
                throw py::index_error();
            m[i] = v;
         })

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                /* Pointer to buffer */
                sizeof(Scalar),          /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::value(),
                1,                       /* Number of dimensions */
                { (size_t) m.size() },   /* Buffer dimensions */
                { sizeof(Scalar) }       /* Strides (in bytes) for each index */
            );
         })

        /* Static initializers */
        .def_static("Zero", [](size_t n) { return Type(Type::Zero(n)); })
        .def_static("Ones", [](size_t n) { return Type(Type::Ones(n)); })
        .def_static("Constant", [](size_t n, Scalar value) { return Type(Type::Constant(n, value)); });
    return vector;
}

/// Creates Python bindings for a dynamic Eigen order-2 tensor (i.e. a matrix)
template <typename Type>
py::class_<Type> bind_eigen_2(py::module &m, const char *name,
                                py::object parent = py::object()) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> matrix(m, name, parent);

    matrix
        /* Constructors */
        .def(py::init<>())
        .def(py::init<size_t, size_t>())
        .def("__init__", [](Type &m, Scalar f) {
            new (&m) Type(1, 1);
            m(0, 0) = f;
        })
        .def("__init__", [](Type &m, py::buffer b) {
            py::buffer_info info = b.request();
            if (info.format != py::format_descriptor<Scalar>::value())
                throw std::runtime_error("Incompatible buffer format!");
            if (info.ndim == 1) {
                new (&m) Type(info.shape[0], 1);
                memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
            } else if (info.ndim == 2) { 
                if (info.strides[0] == sizeof(Scalar)) {
                    new (&m) Type(info.shape[0], info.shape[1]);
                    memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
                } else {
                    new (&m) Type(info.shape[1], info.shape[0]);
                    memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
                    m.transposeInPlace();
                }
            } else {
                throw std::runtime_error("Incompatible buffer dimension!");
            }
        })

        /* Size query functions */
        .def("size", &Type::size)
        .def("cols", &Type::cols)
        .def("rows", &Type::rows)

        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setIdentity", [](Type &m) { m.setIdentity(); })
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })

        /* Resizing */
        .def("resize", [](Type &m, size_t s0, size_t s1) { m.resize(s0, s1); })
        .def("resizeLike", [](Type &m, const Type &m2) { m.resizeLike(m2); })
        .def("conservativeResize", [](Type &m, size_t s0, size_t s1) { m.conservativeResize(s0, s1); })

        /* Component-wise operations */
        .def("cwiseAbs", &Type::cwiseAbs)
        .def("cwiseAbs2", &Type::cwiseAbs2)
        .def("cwiseSqrt", &Type::cwiseSqrt)
        .def("cwiseInverse", &Type::cwiseInverse)
        .def("cwiseMin", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMin(m2); })
        .def("cwiseMax", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMax(m2); })
        .def("cwiseMin", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMin(s); })
        .def("cwiseMax", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMax(s); })
        .def("cwiseProduct", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseProduct(m2); })
        .def("cwiseQuotient", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseQuotient(m2); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * py::self)
        .def_cast(py::self * Scalar())
        .def_cast(py::self / Scalar())

        /* Arithmetic in-place operators */
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)

        .def("transposeInPlace", [](Type &m) { m.transposeInPlace(); })
        /* Other transformations */
        .def("transpose", [](Type &m) -> Type { return m.transpose(); })
        /* Python protocol implementations */
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v;
            return oss.str();
        })
        .def("__getitem__", [](const Type &m, std::pair<size_t, size_t> i) {
            if (i.first >= (size_t) m.rows() || i.second >= (size_t) m.cols())
                throw py::index_error();
            return m(i.first, i.second);
         })
        .def("__setitem__", [](Type &m, std::pair<size_t, size_t> i, Scalar v) {
            if (i.first >= (size_t) m.rows() || i.second >= (size_t) m.cols())
                throw py::index_error();
            m(i.first, i.second) = v;
         })

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                /* Pointer to buffer */
                sizeof(Scalar),          /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::value(),
                2,                       /* Number of dimensions */
                { (size_t) m.rows(),     /* Buffer dimensions */
                  (size_t) m.cols() },
                { sizeof(Scalar),        /* Strides (in bytes) for each index */
                  sizeof(Scalar) * m.rows() }
            );
         })

        /* Static initializers */
        .def_static("Zero", [](size_t n, size_t m) { return Type(Type::Zero(n, m)); })
        .def_static("Ones", [](size_t n, size_t m) { return Type(Type::Ones(n, m)); })
        .def_static("Constant", [](size_t n, size_t m, Scalar value) { return Type(Type::Constant(n, m, value)); })
        .def_static("Identity", [](size_t n, size_t m) { return Type(Type::Identity(n, m)); });
    return matrix;
}

void python_export_vector(py::module &m) {
    bind_eigen_1<VectorX> (m, "VectorX");
    bind_eigen_2<MatrixXd>(m, "MatrixX");

    /* Bindings for <color.h> */
    py::class_<Color3>(m, "Color3")
        .def(py::init<>())
        .def(py::init<Float>())
        .def(py::init<Float, Float, Float>())
        .def("__init__", [](Color3 &c, const std::vector<Float> &v2) {
            if (v2.size() != 3)
                throw std::runtime_error("Incompatible size!");
            memcpy(c.data(), &v2[0], sizeof(Float) * 3);
        })
        .def("__init__", [](Color3 &c, py::buffer b) {
            init_fixed_from_buffer(c, b);
        })
        .def("isValid", &Color3::isValid)
        .def("clamp", &Color3::clamp)
        .def("toLinearRGB", &Color3::toLinearRGB)
        .def("toSRGB", &Color3::toSRGB)
        .def("getLuminance", &Color3::getLuminance)
        .def("__len__", [](const Color3 &) { return 3; })
        .def("__repr__", [](const Color3 &c) {
            std::ostringstream oss;
            oss << c;
            return oss.str();
        })
        .def("__getitem__", [](const Color3 &c, int i) {
            if (i < 0 || i >= 3)
                throw py::index_error();
            return c[i];
         })
        .def("__setitem__", [](Color3 &c, int i, Float v) {
             if (i < 0 || i >= 3)
                 throw py::index_error();
            c[i] = v;
         })
        .def_property("r", [](const Color3 &c) -> Float { return c.r(); },
                           [](Color3 &c, Float v) { c.r() = v; }, "Red channel")
        .def_property("g", [](const Color3 &c) -> Float { return c.g(); },
                           [](Color3 &c, Float v) { c.g() = v; }, "Green channel")
        .def_property("b", [](const Color3 &c) -> Float { return c.b(); },
                           [](Color3 &c, Float v) { c.b() = v; }, "Blue channel");

    /* Bindings for <vector.h> */
    auto vector3 = bind_eigen_1_3<Vector3>(m, "Vector3");
    vector3
        .def("norm", [](const Vector3 &v) { return v.norm(); })
        .def("squaredNorm", [](const Vector3 &v) { return v.squaredNorm(); })
        .def("normalize", [](Vector3 &v) { v.normalize(); })
        .def("normalized", [](const Vector3 &v) -> Vector3 { return v.normalized(); })
        .def("dot", [](const Vector3 &v1, const Vector3 &v2) { return v1.dot(v2); })
        .def("cross", [](const Vector3 &v1, const Vector3 &v2) -> Vector3 { return v1.cross(v2); })
        .def_property("x", [](const Vector3 &v) -> Float { return v.x(); },
                           [](Vector3 &v, Float x) { v.x() = x; }, "X coordinate")
        .def_property("y", [](const Vector3 &v) -> Float { return v.y(); },
                           [](Vector3 &v, Float y) { v.y() = y; }, "Y coordinate")
        .def_property("z", [](const Vector3 &v) -> Float { return v.z(); },
                           [](Vector3 &v, Float z) { v.z() = z; }, "Z coordinate");

    m.def("coordinateSystem", [](const Vector3 &n) {
        Vector3 s, t;
        coordinateSystem(n, s, t);
        return std::make_pair(s, t);
    });

    py::class_<Normal3>(m, "Normal3", vector3)
        /* Constructors */
        .def(py::init<>())
        .def(py::init<Float>())
        .def(py::init<Float, Float, Float>())
        .def("__init__", [](Normal3 &v, const std::vector<Float> &v2) {
            if (v2.size() != Normal3::Dimension)
                throw std::runtime_error("Incompatible size!");
            memcpy(v.data(), &v2[0], sizeof(Float) * Normal3::Dimension);
        })
        .def("__init__", [](Normal3 &v, py::buffer b) {
            init_fixed_from_buffer(v, b);
        });

    bind_eigen_1_3<Point3>(m, "Point3")
        .def_property("x", [](const Point3 &p) -> Float { return p.x(); },
                           [](Point3 &p, Float x) { p.x() = x; }, "X coordinate")
        .def_property("y", [](const Point3 &p) -> Float { return p.y(); },
                           [](Point3 &p, Float y) { p.y() = y; }, "Y coordinate")
        .def_property("z", [](const Point3 &p) -> Float { return p.z(); },
                           [](Point3 &p, Float z) { p.z() = z; }, "Z coordinate");

    py::class_<Frame3>(m, "Frame3", D(TFrame3))
        .def(py::init<>())
        .def(py::init<Frame3>())
        .def(py::init<Vector3>())
        .def(py::init<Vector3, Vector3, Vector3>())
        .def("toLocal", &Frame3::toLocal, D(TFrame3, toLocal))
        .def("toWorld", &Frame3::toWorld, D(TFrame3, toWorld))
        .def("__repr__", [](const Frame3 &f) {
            std::ostringstream oss;
            oss << f;
            return oss.str();
        })
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def_static("sinTheta2", &Frame3::sinTheta2, D(TFrame3, sinTheta2))
        .def_static("cosTheta2", &Frame3::cosTheta2, D(TFrame3, cosTheta2))
        .def_static("tanTheta2", &Frame3::tanTheta2, D(TFrame3, tanTheta2))
        .def_static("sinTheta", &Frame3::sinTheta, D(TFrame3, sinTheta))
        .def_static("cosTheta", &Frame3::cosTheta, D(TFrame3, cosTheta))
        .def_static("tanTheta", &Frame3::tanTheta, D(TFrame3, tanTheta))
        .def_static("sinPhi2", &Frame3::sinPhi2, D(TFrame3, sinPhi2))
        .def_static("cosPhi2", &Frame3::cosPhi2, D(TFrame3, cosPhi2))
        .def_static("sinPhi", &Frame3::sinPhi, D(TFrame3, sinPhi))
        .def_static("cosPhi", &Frame3::cosPhi, D(TFrame3, cosPhi))
        .def_readwrite("s", &Frame3::s, D(TFrame3, s))
        .def_readwrite("t", &Frame3::t, D(TFrame3, t))
        .def_readwrite("n", &Frame3::n, D(TFrame3, n));

    m.attr("Vector") = m.attr("Vector3");
    m.attr("Point")  = m.attr("Point3");
    m.attr("Normal") = m.attr("Normal3");
    m.attr("Frame") = m.attr("Frame3");
    m.attr("Color") = m.attr("Color3");

    py::implicitly_convertible<py::buffer, Normal3>();
    py::implicitly_convertible<py::buffer, Vector3>();
    py::implicitly_convertible<py::buffer, Point3>();
    py::implicitly_convertible<py::buffer, Color3>();
    py::implicitly_convertible<py::buffer, VectorX>();
    py::implicitly_convertible<py::buffer, MatrixX>();

    py::implicitly_convertible<double, Normal3>();
    py::implicitly_convertible<double, Vector3>();
    py::implicitly_convertible<double, Point3>();
    py::implicitly_convertible<double, Color3>();
    py::implicitly_convertible<double, VectorX>();
    py::implicitly_convertible<double, MatrixX>();
}
