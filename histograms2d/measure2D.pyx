# distutils: language = c++

from DOT_Measure2D cimport Measure2D as PyMeasure2D


cdef class Measure2D:
    cdef PyMeasure2D mu

    def __cinit__(self):
        self.mu = PyMeasure2D()

    def add(self, i, j, w):
        return self.mu.add(i, w)

    def size(self):
        return self.mu.size()