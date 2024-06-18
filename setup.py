from setuptools.extension import Extension
from setuptools import setup
from Cython.Build import cythonize

setup(
  name = "rayx",
  ext_modules = cythonize([
    Extension("rayx", 
      sources=["rayx.pyx"],
      language="c++",
      include_dirs=[
        "rayx/Extern/glm/glm",
        "rayx/Extern/rapidxml-1.13",
        "rayx/Intern/rayx-core/src",
      ],
      libraries=["rayx-core"],
      library_dirs=["rayx/build/Intern/rayx-core"],
      extra_compile_args=["-std=c++20"],
      extra_link_args=["-std=c++20"],
    )
  ])
)