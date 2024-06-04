# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp.string cimport string
from libcpp.vector cimport vector

import pandas as pd

cdef extern from "HelloCpp.cpp":
  pass


cdef extern from "HelloCpp.h":
  cdef cppclass HelloCpp:
    HelloCpp() except +
    HelloCpp(string) except +
    void sayHello()
    void setRecipient(string)

  cdef cppclass VectorSource:
    VectorSource() except +
    VectorSource(int) except +
    void set_size(int)
    vector[int] getData()

cdef class Hello:
  cdef HelloCpp c_hello
  def __cinit__(self, str name):
    self.c_hello = HelloCpp(name)

  def say_hello(self):
    self.c_hello.sayHello()

  def set_recipient(self, str name):
    self.c_hello.setRecipient(name)

cdef class Vector:
  cdef VectorSource c_vector
  def __cinit__(self, int size):
    self.c_vector = VectorSource(size)

  def set_size(self, int size):
    self.c_vector.set_size(size)

  def get_data(self):
    cdef vector[int] data = self.c_vector.getData()
    return pd.DataFrame(data, columns=['data'])