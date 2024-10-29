from setuptools import setup, Extension, find_namespace_packages
from Cython.Build import cythonize
import os

# download boost
if not os.path.exists("boost_1_86_0"):
  os.system("curl -LO https://boostorg.jfrog.io/artifactory/main/release/1.86.0/source/boost_1_86_0.tar.gz")
  os.system("tar -xzf boost_1_86_0.tar.gz")

boost_dir = os.path.abspath("boost_1_86_0")

os.system(f"mkdir -p rayx/build && cd rayx/build && cmake -DRAYX_STATIC_LIB=ON -DBoost_INCLUDE_DIR={boost_dir} .. && make -j")
os.system("mkdir -p pkgs/rayx-data/share/RAYX")
os.system("cp -r rayx/Data pkgs/rayx-data/share/RAYX")

extensions = [
  Extension("rayx", 
    ["rayx.pyx"],
    language="c++",
    include_dirs=[
      "rayx/Extern/glm/glm",
      "rayx/Extern/rapidxml-1.13",
      "rayx/Extern/alpaka/include",
      "rayx/Intern/rayx-core/src",
      "boost_1_86_0"
    ],
    libraries=["gomp"],
    extra_compile_args=["-std=c++20", "-fopenmp"],
    extra_link_args=["-std=c++20", "-fopenmp", "-static-libstdc++"],
    extra_objects=["rayx/build/Intern/rayx-core/librayx-core.a"],
  )
]

setup(
  name = "rayx",
  ext_modules = cythonize(extensions),
  packages=find_namespace_packages(where="pkgs"),
  package_dir={"": "pkgs"},
  package_data={
    "rayx-data.Data.nff": ["*.nff"],
    "rayx-data.Data.PALIK": ["*.NKP"],
  }
)