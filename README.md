# ``layerlab``: A computational toolbox for layered materials

``layerlab`` is a Python-based toolbox for computations
involving layered materials that implements a model described in the paper
[A Comprehensive Framework for Rendering Layered Materials](http://www.cs.cornell.edu/projects/layered-sg14/)
by Wenzel Jakob, Eugene D'Eon, Otto Jakob and Steve Marschner.

A layered BSDF model describes the directional reflectance properties of a
material whose internal structure consists of a stack of scattering and/or
absorbing layers separated by smooth or rough interfaces. The bottom of the
stack could be an opaque interface (such as a metal) or a transparent one. Such
structural decompositions into layers and interfaces dramatically enlarge the
size of the “language” that is available to describe materials, and for this
reason they have been the focus of considerable interest in computer graphics
in the last years.

See [http://www.mitsuba-renderer.org/~wenzel/papers/layerlab-sg15.pdf](http://www.mitsuba-renderer.org/~wenzel/papers/layerlab-sg15.pdf)
for a tutorial on using ``layerlab``.

A few recipe-style examples are available in the ``recipes`` directory.
