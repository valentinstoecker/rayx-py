# cython: c_string_type=unicode, c_string_encoding=utf8

import pandas as pd

from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp.memory cimport shared_ptr
from cython.operator cimport dereference

cdef extern from "rayx/Extern/glm/glm/glm.hpp" namespace "glm":
    cdef cppclass dvec3:
        dvec3() except +
        dvec3(double x, double y, double z) except +
        double x, y, z
        double& operator[](int)
    
    cdef cppclass dvec4:
        dvec4() except +
        dvec4(double x, double y, double z, double w) except +
        double x, y, z, w
        double& operator[](int)
    
    cdef cppclass dmat4x4:
        dmat4x4() except +
        dvec4& operator[](int)

cdef extern from "rayx/Intern/rayx-core/raycorePCH.h":
    pass

cdef extern from "rayx/Intern/rayx-core/src/Shader/Ray.h" namespace "RAYX":
    cdef cppclass Ray:
        Ray() except +
        dvec3 m_position

        double m_eventType

        dvec3 m_direction

        double m_energy

        dvec4 m_stokes

        double m_pathLength

        double m_order

        double m_lastElement

        double m_sourceID

cdef extern from "rayx/Intern/rayx-core/src/angle.h":
    cdef cppclass Rad:
        Rad() except +
        Rad(double) except +
        double rad

cdef extern from "rayx/Intern/rayx-core/src/Data/xml.h" namespace "RAYX":
    cdef cppclass Misalignment:
        Misalignment() except +
        double m_translationXerror, m_translationYerror, m_translationZerror
        Rad m_rotationXerror, m_rotationYerror, m_rotationZerror
cdef extern from "rayx/Intern/rayx-core/src/Shader/Element.h" namespace "RAYX":
    cdef cppclass SlopeError:
        SlopeError() except +
        double m_sag, m_mer, m_thermalDistortionAmp, m_thermalDistortionSigmaX, m_thermalDistortionSigmaZ, m_cylindricalBowingAmp, m_cylindricalBowingRadius

cdef extern from "rayx/Intern/rayx-core/src/Material/Material.h" namespace "RAYX":
    cpdef enum class Material(int):
        REFLECTIVE = -2,
        VACUUM = -1,
        H = 1,
        He = 2,
        Li = 3,
        Be = 4,
        B = 5,
        C = 6,
        N = 7,
        O = 8,
        F = 9,
        Ne = 10,
        Na = 11,
        Mg = 12,
        Al = 13,
        Si = 14,
        P = 15,
        S = 16,
        Cl = 17,
        Ar = 18,
        K = 19,
        Ca = 20,
        Sc = 21,
        Ti = 22,
        V = 23,
        Cr = 24,
        Mn = 25,
        Fe = 26,
        Co = 27,
        Ni = 28,
        Cu = 29,
        Zn = 30,
        Ga = 31,
        Ge = 32,
        As = 33,
        Se = 34,
        Br = 35,
        Kr = 36,
        Rb = 37,
        Sr = 38,
        Y = 39,
        Zr = 40,
        Nb = 41,
        Mo = 42,
        Tc = 43,
        Ru = 44,
        Rh = 45,
        Pd = 46,
        Ag = 47,
        Cd = 48,
        In = 49,
        Sn = 50,
        Sb = 51,
        Te = 52,
        I = 53,
        Xe = 54,
        Cs = 55,
        Ba = 56,
        La = 57,
        Ce = 58,
        Pr = 59,
        Nd = 60,
        Pm = 61,
        Sm = 62,
        Eu = 63,
        Gd = 64,
        Tb = 65,
        Dy = 66,
        Ho = 67,
        Er = 68,
        Tm = 69,
        Yb = 70,
        Lu = 71,
        Hf = 72,
        Ta = 73,
        W = 74,
        Re = 75,
        Os = 76,
        Ir = 77,
        Pt = 78,
        Au = 79,
        Hg = 80,
        Tl = 81,
        Pb = 82,
        Bi = 83,
        Po = 84,
        At = 85,
        Rn = 86,
        Fr = 87,
        Ra = 88,
        Ac = 89,
        Th = 90,
        Pa = 91,
        U = 92,

cdef extern from "rayx/Intern/rayx-core/src/Shader/Cutout.h":
    cdef cppclass Cutout:
        Cutout() except +
        double m_type
    
    cdef cppclass RectCutout:
        RectCutout() except +
        double m_width, m_length
    
    Cutout serializeRect(RectCutout)
    RectCutout deserializeRect(Cutout)

    cdef cppclass EllipticalCutout:
        EllipticalCutout() except +
        double m_diameter_x, m_diameter_z
    
    Cutout serializeElliptical(EllipticalCutout)
    EllipticalCutout deserializeElliptical(Cutout)

    cdef cppclass TrapezoidCutout:
        TrapezoidCutout() except +
        double m_widthA, m_widthB, m_length

    Cutout serializeTrapezoid(TrapezoidCutout)
    TrapezoidCutout deserializeTrapezoid(Cutout)

    Cutout serializeUnlimited()

cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass array6d "std::array<double, 6>":
        array6d() except +
        double& operator[](int)

cdef extern from "rayx/Intern/rayx-core/src/Shader/Surface.h":
    cdef cppclass Surface:
        Surface() except +
        double m_type
    
    cdef cppclass QuadricSurface:
        QuadricSurface() except +
        int m_icurv
        double m_a11, m_a12, m_a13, m_a14, m_a22, m_a23, m_a24, m_a33, m_a34, m_a44
    
    Surface serializeQuadric(QuadricSurface)
    QuadricSurface deserializeQuadric(Surface)

    cdef cppclass ToroidSurface:
        ToroidSurface() except +
        double m_longRadius, m_shortRadius, m_toroidType
    
    Surface serializeToroid(ToroidSurface)
    ToroidSurface deserializeToroid(Surface)

    Surface serializePlaneXZ()

    cdef cppclass CubicSurface:
        CubicSurface() except +
        int m_icurv
        double m_a11, m_a12, m_a13, m_a14, m_a22, m_a23, m_a24, m_a33, m_a34, m_a44, m_b12, m_b13, m_b21, m_b23, m_b31, m_b32, m_psi

    Surface serializeCubic(CubicSurface)
    CubicSurface deserializeCubic(Surface)

cdef extern from "rayx/Intern/rayx-core/src/Beamline/Objects/BehaviourType.h" namespace "RAYX":
    cpdef enum class CentralBeamstop(int):
        None,
        Rectangle,
        Elliptical,

cdef extern from "rayx/Intern/rayx-core/src/Beamline/Objects/SurfaceType.h" namespace "RAYX":
    cpdef enum class CylinderDirection(int):
        LongRadiusR, 
        ShortRadiusRho,

cdef extern from "rayx/Intern/rayx-core/src/Beamline/OpticalElement.h" namespace "RAYX":
    cpdef enum class FigureRotation(int):
        Yes,
        Plane,
        A11,

    cpdef enum class CurvatureType(int):
        Plane,
        Toroidal,
        Spherical,
        Cubic,
        Cone,
        Cylinder,
        Ellipsoid,
        Paraboloid,
        Quadric,
        RzpSphere,
    
    cpdef enum class BehaviourType(int):
        Mirror,
        Grating,
        Slit,
        Rzp,
        ImagePlane,

cdef extern from "rayx/Intern/rayx-core/src/Beamline/LightSource.h" namespace "RAYX":
    cpdef enum class SourceDist(int):
        Uniform,
        Gaussian,
        Thirds,
        Circle,
    
    cpdef enum class SpreadType(int):
        HardEdge,
        SoftEdge,
        SeperateEnergies,

    cpdef enum class EnergyDistributionType(int):
        File,
        Values,
        Total,
        Param,

    cpdef enum class EnergySpreadUnit(int):
        EU_PERCENT,
        EU_eV,

    cpdef enum class ElectronEnergyOrientation(int):
        Clockwise,
        Counterclockwise,

    cpdef enum class SigmaType(int):
        ST_STANDARD,
        ST_ACCURATE,

cdef extern from "rayx/Intern/rayx-core/src/DesignElement/Value.h" namespace "RAYX":
    cpdef enum class ValueType(int):
        Undefined,
        Double,
        Int,
        String,
        Map,
        Dvec4,
        Dmat4x4,
        Bool,
        Rad,
        Material,
        Misalignment,
        CentralBeamStop,
        Cutout,
        CylinderDirection,
        FigureRotation,
        CurvatureType,
        Surface,
        SourceDist,
        SpreadType,
        EnergyDistributionType,
        EnergySpreadUnit,
        ElectronEnergyOrientation,
        SigmaType,
        BehaviourType,
    
    cdef cppclass DesignMap:
        DesignMap() except +

        void operator=(double)
        void operator=(int)
        void operator=(bool)
        void operator=(string)
        void operator=(unordered_map[string, shared_ptr[DesignMap]])
        void operator=(dvec4)
        void operator=(dmat4x4)
        void operator=(Rad)
        void operator=(Material)
        void operator=(Misalignment)
        void operator=(CentralBeamstop)
        void operator=(Cutout)
        void operator=(CylinderDirection)
        void operator=(FigureRotation)
        void operator=(CurvatureType)
        void operator=(Surface)
        void operator=(SourceDist)
        void operator=(SpreadType)
        void operator=(EnergyDistributionType)
        void operator=(EnergySpreadUnit)
        void operator=(ElectronEnergyOrientation)
        void operator=(SigmaType)
        void operator=(BehaviourType)

        ValueType type()

        double as_double()
        int as_int()
        cpp_bool as_bool()
        string as_string()
        unordered_map[string, shared_ptr[DesignMap]] as_map()
        dvec4 as_dvec4()
        dmat4x4 as_dmat4x4()
        Rad as_rad()
        Material as_material()
        Misalignment as_misalignment()
        CentralBeamstop as_centralBeamStop()
        Cutout as_cutout()
        CylinderDirection as_cylinderDirection()
        FigureRotation as_figureRotation()
        CurvatureType as_curvatureType()
        Surface as_surface()
        SourceDist as_sourceDist()
        SpreadType as_energySpreadType()
        EnergyDistributionType as_energyDistType()
        EnergySpreadUnit as_energySpreadUnit()
        ElectronEnergyOrientation as_electronEnergyOrientation()
        SigmaType as_sigmaType()
        BehaviourType as_behaviourType()

        DesignMap& operator[](string)


cdef extern from "rayx/Intern/rayx-core/src/DesignElement/DesignElement.h" namespace "RAYX":
    cdef cppclass DesignElement:
        DesignElement() except +

        DesignMap m_elementParameters

cdef extern from "rayx/Intern/rayx-core/src/DesignElement/DesignSource.h" namespace "RAYX":
    cdef cppclass DesignSource:
        DesignSource() except +

        DesignMap m_elementParameters

cdef extern from "rayx/Intern/rayx-core/src/Beamline/Beamline.h" namespace "RAYX":
    cdef cppclass Beamline:
        Beamline() except +

        vector[DesignElement] m_DesignElements
        vector[DesignSource] m_DesignSources

cdef extern from "rayx/Intern/rayx-core/src/Data/Importer.h" namespace "RAYX":
    cdef Beamline importBeamline(string path) except +

cdef extern from "rayx/Intern/rayx-core/src/Tracer/CpuTracer.h" namespace "RAYX":
    cpdef enum class Sequential(int):
        No,
        Yes,
    
    cdef cppclass CpuTracer:
        CpuTracer() except +
        vector[vector[Ray]] trace(Beamline&, Sequential, int, int, int)

#cdef extern from "rayx/Intern/rayx-core/src/Tracer/VulkanTracer.h" namespace "RAYX":
#    cdef cppclass VulkanTracer:
#        VulkanTracer() except +
#        vector[vector[Ray]] trace(Beamline&, Sequential, int, int, int)

cdef designMapToPy(DesignMap dm):
    cdef dvec4 dv
    cdef dmat4x4 dm4
    cdef Misalignment mis
    cdef RectCutout rect
    cdef EllipticalCutout ell
    cdef TrapezoidCutout trap
    cdef QuadricSurface quad
    cdef ToroidSurface tor
    cdef CubicSurface cub
    if dm.type() == ValueType.Double:
        return dm.as_double()
    elif dm.type() == ValueType.Int:
        return dm.as_int()
    elif dm.type() == ValueType.Bool:
        return dm.as_bool()
    elif dm.type() == ValueType.String:
        return dm.as_string()
    elif dm.type() == ValueType.Map:
        py_dm: dict = {}
        for p in dm.as_map():
            py_dm[p.first] = designMapToPy(dereference(p.second))
        return py_dm
    elif dm.type() == ValueType.Dvec4:
        dv = dm.as_dvec4()
        return [dv[0], dv[1], dv[2], dv[3]]
    elif dm.type() == ValueType.Dmat4x4:
        dm4 = dm.as_dmat4x4()
        return [[dm4[0][0], dm4[0][1], dm4[0][2], dm4[0][3]],
                [dm4[1][0], dm4[1][1], dm4[1][2], dm4[1][3]],
                [dm4[2][0], dm4[2][1], dm4[2][2], dm4[2][3]],
                [dm4[3][0], dm4[3][1], dm4[3][2], dm4[3][3]]]
    elif dm.type() == ValueType.Rad:
        return {
            "_type": "Rad",
            "value": dm.as_rad().rad
        }
    elif dm.type() == ValueType.Material:
        return dm.as_material()
    elif dm.type() == ValueType.Misalignment:
        mis = dm.as_misalignment()
        return {
            "_type": "Misalignment",
            "translationXerror": mis.m_translationXerror,
            "translationYerror": mis.m_translationYerror,
            "translationZerror": mis.m_translationZerror,
            "rotationXerror": mis.m_rotationXerror.rad,
            "rotationYerror": mis.m_rotationYerror.rad,
            "rotationZerror": mis.m_rotationZerror.rad
        }
    elif dm.type() == ValueType.CentralBeamStop:
        return dm.as_centralBeamStop()
    elif dm.type() == ValueType.Cutout:
        if dm.as_cutout().m_type == 0:
            rect = deserializeRect(dm.as_cutout())
            return {
                "_type": "RectCutout",
                "width": rect.m_width,
                "length": rect.m_length
            }
        elif dm.as_cutout().m_type == 1:
            ell = deserializeElliptical(dm.as_cutout())
            return {
                "_type": "EllipticalCutout",
                "diameter_x": ell.m_diameter_x,
                "diameter_z": ell.m_diameter_z
            }
        elif dm.as_cutout().m_type == 2:
            trap = deserializeTrapezoid(dm.as_cutout())
            return {
                "_type": "TrapezoidCutout",
                "widthA": trap.m_widthA,
                "widthB": trap.m_widthB,
                "length": trap.m_length
            }
        elif dm.as_cutout().m_type == 3:
            return {
                "_type": "UnlimitedCutout"
            }
    elif dm.type() == ValueType.CylinderDirection:
        return dm.as_cylinderDirection()
    elif dm.type() == ValueType.FigureRotation:
        return dm.as_figureRotation()
    elif dm.type() == ValueType.CurvatureType:
        return dm.as_curvatureType()
    elif dm.type() == ValueType.Surface:
        if dm.as_surface().m_type == 0:
            quad = deserializeQuadric(dm.as_surface())
            return {
                "_type": "QuadricSurface",
                "icurv": quad.m_icurv,
                "a11": quad.m_a11,
                "a12": quad.m_a12,
                "a13": quad.m_a13,
                "a14": quad.m_a14,
                "a22": quad.m_a22,
                "a23": quad.m_a23,
                "a24": quad.m_a24,
                "a33": quad.m_a33,
                "a34": quad.m_a34,
                "a44": quad.m_a44
            }
        elif dm.as_surface().m_type == 1:
            tor = deserializeToroid(dm.as_surface())
            return {
                "_type": "ToroidSurface",
                "longRadius": tor.m_longRadius,
                "shortRadius": tor.m_shortRadius,
                "toroidType": tor.m_toroidType
            }
        elif dm.as_surface().m_type == 2:
            return {
                "_type": "PlaneXZSurface"
            }
        elif dm.as_surface().m_type == 3:
            cub = deserializeCubic(dm.as_surface())
            return {
                "_type": "CubicSurface",
                "icurv": cub.m_icurv,
                "a11": cub.m_a11,
                "a12": cub.m_a12,
                "a13": cub.m_a13,
                "a14": cub.m_a14,
                "a22": cub.m_a22,
                "a23": cub.m_a23,
                "a24": cub.m_a24,
                "a33": cub.m_a33,
                "a34": cub.m_a34,
                "a44": cub.m_a44,
                "b12": cub.m_b12,
                "b13": cub.m_b13,
                "b21": cub.m_b21,
                "b23": cub.m_b23,
                "b31": cub.m_b31,
                "b32": cub.m_b32,
                "psi": cub.m_psi
            }
    elif dm.type() == ValueType.SourceDist:
        return dm.as_sourceDist()
    elif dm.type() == ValueType.SpreadType:
        return dm.as_energySpreadType()
    elif dm.type() == ValueType.EnergyDistributionType:
        return dm.as_energyDistType()
    elif dm.type() == ValueType.EnergySpreadUnit:
        return dm.as_energySpreadUnit()
    elif dm.type() == ValueType.ElectronEnergyOrientation:
        return dm.as_electronEnergyOrientation()
    elif dm.type() == ValueType.SigmaType:
        return dm.as_sigmaType()
    elif dm.type() == ValueType.BehaviourType:
        return dm.as_behaviourType()
    else:
        return None

cdef DesignMap dictToDesignMap(d):
    cdef DesignMap dm = DesignMap()
    cdef DesignMap* inner_dm
    cdef unordered_map[string, shared_ptr[DesignMap]] map
    cdef Misalignment mis

    cdef Cutout cut
    cdef RectCutout rect
    cdef EllipticalCutout ell
    cdef TrapezoidCutout trap

    cdef Surface surf
    cdef QuadricSurface quad
    cdef ToroidSurface tor
    cdef CubicSurface cub


    cdef double db
    cdef int i
    cdef cpp_bool b
    cdef string cstring

    cdef Material mat
    cdef CentralBeamstop cb
    cdef CylinderDirection cd
    cdef FigureRotation fr
    cdef CurvatureType ct
    cdef SourceDist sd
    cdef SpreadType st
    cdef EnergyDistributionType ed
    cdef EnergySpreadUnit es
    cdef ElectronEnergyOrientation eo
    cdef SigmaType sig
    cdef BehaviourType bt

    if isinstance(d, dict):
        dtype = d.get("_type")
        if dtype == None:
            for k, v in d.items():
                inner_dm = new DesignMap()
                inner_dm[0] = dictToDesignMap(v)
                map[k] = shared_ptr[DesignMap](inner_dm)
            dm = map
        elif dtype == "Rad":
            dm = Rad(d.get("value"))
        elif dtype == "Misalignment":
            mis = Misalignment()
            mis.m_translationXerror = d.get("translationXerror")
            mis.m_translationYerror = d.get("translationYerror")
            mis.m_translationZerror = d.get("translationZerror")
            mis.m_rotationXerror = Rad(d.get("rotationXerror"))
            mis.m_rotationYerror = Rad(d.get("rotationYerror"))
            mis.m_rotationZerror = Rad(d.get("rotationZerror"))
            dm = mis
        elif dtype == "RectCutout":
            rect = RectCutout()
            rect.m_width = d.get("width")
            rect.m_length = d.get("length")
            cut = serializeRect(rect)
            dm = cut
        elif dtype == "EllipticalCutout":
            ell = EllipticalCutout()
            ell.m_diameter_x = d.get("diameter_x")
            ell.m_diameter_z = d.get("diameter_z")
            cut = serializeElliptical(ell)
            dm = cut
        elif dtype == "TrapezoidCutout":
            trap = TrapezoidCutout()
            trap.m_widthA = d.get("widthA")
            trap.m_widthB = d.get("widthB")
            trap.m_length = d.get("length")
            cut = serializeTrapezoid(trap)
            dm = cut
        elif dtype == "UnlimitedCutout":
            dm = serializeUnlimited()
        elif dtype == "QuadricSurface":
            quad = QuadricSurface()
            quad.m_icurv = d.get("icurv")
            quad.m_a11 = d.get("a11")
            quad.m_a12 = d.get("a12")
            quad.m_a13 = d.get("a13")
            quad.m_a14 = d.get("a14")
            quad.m_a22 = d.get("a22")
            quad.m_a23 = d.get("a23")
            quad.m_a24 = d.get("a24")
            quad.m_a33 = d.get("a33")
            quad.m_a34 = d.get("a34")
            quad.m_a44 = d.get("a44")
            surf = serializeQuadric(quad)
            dm = surf
        elif dtype == "ToroidSurface":
            tor = ToroidSurface()
            tor.m_longRadius = d.get("longRadius")
            tor.m_shortRadius = d.get("shortRadius")
            tor.m_toroidType = d.get("toroidType")
            surf = serializeToroid(tor)
            dm = surf
        elif dtype == "PlaneXZSurface":
            surf = serializePlaneXZ()
            dm = surf
        elif dtype == "CubicSurface":
            cub = CubicSurface()
            cub.m_icurv = d.get("icurv")
            cub.m_a11 = d.get("a11")
            cub.m_a12 = d.get("a12")
            cub.m_a13 = d.get("a13")
            cub.m_a14 = d.get("a14")
            cub.m_a22 = d.get("a22")
            cub.m_a23 = d.get("a23")
            cub.m_a24 = d.get("a24")
            cub.m_a33 = d.get("a33")
            cub.m_a34 = d.get("a34")
            cub.m_a44 = d.get("a44")
            cub.m_b12 = d.get("b12")
            cub.m_b13 = d.get("b13")
            cub.m_b21 = d.get("b21")
            cub.m_b23 = d.get("b23")
            cub.m_b31 = d.get("b31")
            cub.m_b32 = d.get("b32")
            cub.m_psi = d.get("psi")
            surf = serializeCubic(cub)
            dm = surf
    elif isinstance(d, list):
        if len(d) == 4:
            if isinstance(d[0], list):
                dm4 = dmat4x4()
                for i in range(4):
                    for j in range(4):
                        dm4[i][j] = d[i][j]
                dm = dm4
            elif isinstance(d[0], float):
                dv = dvec4()
                for i in range(4):
                    dv[i] = d[i]
                dm = dv
    elif isinstance(d, Material):
        mat = d
        dm = mat
    elif isinstance(d, CentralBeamstop):
        cb = d
        dm = cb
    elif isinstance(d, CylinderDirection):
        cd = d
        dm = cd
    elif isinstance(d, FigureRotation):
        fr = d
        dm = fr
    elif isinstance(d, CurvatureType):
        ct = d
        dm = ct
    elif isinstance(d, SourceDist):
        sd = d
        dm = sd
    elif isinstance(d, SpreadType):
        st = d
        dm = st
    elif isinstance(d, EnergyDistributionType):
        ed = d
        dm = ed
    elif isinstance(d, EnergySpreadUnit):
        es = d
        dm = es
    elif isinstance(d, ElectronEnergyOrientation):
        eo = d
        dm = eo
    elif isinstance(d, SigmaType):
        sig = d
        dm = sig
    elif isinstance(d, BehaviourType):
        bt = d
        dm = bt
    elif isinstance(d, int):
        i = d
        dm = i 
    elif isinstance(d, float):
        db = d
        dm = db
    elif isinstance(d, bool):
        b = d
        dm = b
    elif isinstance(d, str):
        cstring = d
        dm = cstring
    return dm

cdef class BeamlineObj:
    cdef Beamline c_beamline
    def __cinit__(self):
        self.c_beamline = Beamline()
    
    def __init__(self, str path):
        self.c_beamline = importBeamline(path)

    def trace(self, gpu = False, sequential = True, max_batch_size = 100000, thread_count = 0, max_events = -1, start_event_id = 0):
        cdef Sequential seq
        if sequential:
            seq = Sequential.Yes
        else:
            seq = Sequential.No
        cdef int maxEvents = self.c_beamline.m_DesignElements.size() + 2
        if max_events >= 0:
            maxEvents = max_events
        cdef vector[vector[Ray]] rays
        cdef CpuTracer cpu_tracer = CpuTracer()
        #cdef VulkanTracer gpu_tracer = VulkanTracer()
        #if gpu:
            #rays = gpu_tracer.trace(
            #    self.c_beamline,
            #    seq,
            #    max_batch_size, 
            #    thread_count, 
            #    maxEvents
            #)
        #else:
        rays = cpu_tracer.trace(
            self.c_beamline,
            seq,
            max_batch_size, 
            thread_count, 
            maxEvents
        )
        des = self.getDesignElements()
        dss = self.getDesignSources()
        py_rays = []
        for i in range(rays.size()):
            for j in range(rays[i].size()):
                r_obj = RayObj()
                r_obj.c_ray = rays[i][j]
                r_dict = r_obj.to_dict()
                r_dict["ray_index"] = i
                r_dict["event_index"] = j
                r_dict["lastElement"] = des[int(r_dict["lastElement"])]["name"]
                r_dict["sourceID"] = dss[int(r_dict["sourceID"])]["name"]

                py_rays.append(r_dict)
        return pd.DataFrame(py_rays)

    def getDesignElements(self):
        py_elements = []
        for e in self.c_beamline.m_DesignElements:
            py_elements.append(designMapToPy(e.m_elementParameters))
        return py_elements
    
    def setDesignElements(self, elements):
        self.c_beamline.m_DesignElements.clear()
        cdef DesignElement de = DesignElement()
        for e in elements:
            de.m_elementParameters = dictToDesignMap(e)
            self.c_beamline.m_DesignElements.push_back(de)

    def getDesignSources(self):
        py_sources = []
        for s in self.c_beamline.m_DesignSources:
            py_sources.append(designMapToPy(s.m_elementParameters))
        return py_sources
    
    def setDesignSources(self, sources):
        self.c_beamline.m_DesignSources.clear()
        cdef DesignSource ds = DesignSource()
        for s in sources:
            ds.m_elementParameters = dictToDesignMap(s)
            self.c_beamline.m_DesignSources.push_back(ds)

cdef class RayObj:
    cdef Ray c_ray
    def __cinit__(self):
        self.c_ray = Ray()

    def to_dict(self):
        return {
            "position": self.get_position(),
            "eventType": self.get_eventType(),
            "direction": self.get_direction(),
            "energy": self.get_energy(),
            "stokes": self.get_stokes(),
            "pathlength": self.get_pathlength(),
            "order": self.get_order(),
            "lastElement": self.get_lastElement(),
            "sourceID": self.get_sourceID()
        }

    def get_position(self):
        return [self.c_ray.m_position.x, self.c_ray.m_position.y, self.c_ray.m_position.z]
    
    def set_position(self, position):
        self.c_ray.m_position = dvec3(position[0], position[1], position[2])
    
    def get_eventType(self):
        return self.c_ray.m_eventType
    
    def set_eventType(self, eventType):
        self.c_ray.m_eventType = eventType
    
    def get_direction(self):
        return [self.c_ray.m_direction.x, self.c_ray.m_direction.y, self.c_ray.m_direction.z]

    def set_direction(self, direction):
        self.c_ray.m_direction = dvec3(direction[0], direction[1], direction[2])
    
    def get_energy(self):
        return self.c_ray.m_energy
    
    def set_energy(self, energy):
        self.c_ray.m_energy = energy
    
    def get_stokes(self):
        return [self.c_ray.m_stokes.x, self.c_ray.m_stokes.y, self.c_ray.m_stokes.z, self.c_ray.m_stokes.w]
    
    def set_stokes(self, stokes):
        self.c_ray.m_stokes = dvec4(stokes[0], stokes[1], stokes[2], stokes[3])
    
    def get_pathlength(self):
        return self.c_ray.m_pathLength
    
    def set_pathlength(self, pathLength):
        self.c_ray.m_pathLength = pathLength
    
    def get_order(self):
        return self.c_ray.m_order
    
    def set_order(self, order):
        self.c_ray.m_order = order
    
    def get_lastElement(self):
        return self.c_ray.m_lastElement
    
    def set_lastElement(self, lastElement):
        self.c_ray.m_lastElement = lastElement
    
    def get_sourceID(self):
        return self.c_ray.m_sourceID
    
    def set_sourceID(self, sourceID):
        self.c_ray.m_sourceID = sourceID

    def __str__(self):
        s = "Ray:\n"
        s += "Position: " + str(self.get_position()) + "\n"
        s += "Event Type: " + str(self.get_eventType()) + "\n"
        s += "Direction: " + str(self.get_direction()) + "\n"
        s += "Energy: " + str(self.get_energy()) + "\n"
        s += "Stokes: " + str(self.get_stokes()) + "\n"
        s += "Pathlength: " + str(self.get_pathlength()) + "\n"
        s += "Order: " + str(self.get_order()) + "\n"
        s += "Last Element: " + str(self.get_lastElement()) + "\n"
        s += "Source ID: " + str(self.get_sourceID()) + "\n"
        return s
    
    def __repr__(self):
        return self.__str__()
    