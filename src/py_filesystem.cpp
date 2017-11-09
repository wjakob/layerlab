#include <layer/common.h>
#include <filesystem/path.h>
#include <filesystem/resolver.h>
#include <layer/mmap.h>

#include "python.h"

void python_export_filesystem(py::module &m) {
    using namespace fs;

    py::module fs_module =
        m.def_submodule("filesystem", "Cross-platform filesystem support");

    py::class_<path> path_class(fs_module, "path");

    path_class
        .def(py::init<>())
        .def(py::init<const path &>())
        .def(py::init<const std::string &>())
        .def("__len__", &path::length)
        .def("file_size", &path::file_size)
        .def("empty", &path::empty)
        .def("is_absolute", &path::is_absolute)
        .def("make_absolute", &path::make_absolute)
        .def("exists", &path::exists)
        .def("is_directory", &path::is_directory)
        .def("is_file", &path::is_file)
        .def("extension", &path::extension)
        .def("parent_path", &path::parent_path)
        .def("remove_file", &path::remove_file)
        .def("resize_file", &path::resize_file)
        .def("str", [](const path &p) { return p.str(); })
        .def("str", [](const path &p, path::path_type t) { return p.str(t); })
        .def("set", [](path &p, const std::string &v, path::path_type t) { p.set(v, t); })
        .def("set", [](path &p, const std::string &v) { p.set(v); })
        .def(py::self / py::self)
        .def("__repr__", [](const path &p) { return p.str(); })
        .def_static("getcwd", &path::getcwd);

    py::enum_<path::path_type>(path_class, "path_type")
        .value("windows_path", path::windows_path)
        .value("posix_path", path::posix_path)
        .value("native_path", path::native_path)
        .export_values();

    py::class_<resolver>(fs_module, "resolver")
        .def(py::init<>())
        .def("__len__", &resolver::size)
        .def("append", &resolver::append)
        .def("prepend", &resolver::prepend)
        .def("resolve", &resolver::resolve)
        .def("__getitem__", [](const resolver &r, size_t i) {
             if (i >= r.size())
                 throw py::index_error();
             return r[i];
         })
        .def("__setitem__", [](resolver &r, size_t i, path &v) {
             if (i >= r.size())
                 throw py::index_error();
             r[i] = v;
         })
        .def("__delitem__", [](resolver &r, size_t i) {
             if (i >= r.size())
                 throw py::index_error();
             r.erase(r.begin() + i);
         })
        .def("__repr__", [](const resolver &r) {
            std::ostringstream oss;
            oss << r;
            return oss.str();
        });

    py::class_<MemoryMappedFile>(fs_module, "MemoryMappedFile", py::buffer_protocol())
        .def(py::init<fs::path, size_t>())
        .def(py::init<fs::path, bool>())
        .def(py::init<fs::path>())
        .def("resize", &MemoryMappedFile::resize)
        .def("__repr__", &MemoryMappedFile::toString)
        .def("size", &MemoryMappedFile::size)
        .def("filename", &MemoryMappedFile::filename)
        .def("readOnly", &MemoryMappedFile::readOnly)
        .def_static("createTemporary", &MemoryMappedFile::createTemporary)
        .def_buffer([](MemoryMappedFile &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),
                sizeof(uint8_t),
                py::format_descriptor<uint8_t>::format(),
                1,
                { (size_t) m.size() },
                { sizeof(uint8_t) }
            );
        });

    py::implicitly_convertible<std::string, path>();
}
