A little single file library for working with SphericalHarmonicsL2

The interesting functions are SHMatrix, RawSphericalHarmonicsL2 and SHUtility.

# SHTools
A little single file library for working with SphericalHarmonicsL2. The file introduces 3 interesting types:
- `SHUtility`: A static class contain helper functions for working with Spherical Harmonics.
- `RawSphericalHarmonicsL2`: A type that can be implicitly converted to and from SphericalHarmonicsL2, but doing so accounts for the fact that Unity stores SH data in a strange, pre-multiplied convention, which isn't easy to manipulate. Several of the functions in the library work on this type instead of SphericalHarmonicsL2.
- `SHMatrix`: A block matrix representing a rotation of a set of Spherical Harmonics coefficients.

# Function overview

```cs
static class SHUtility {
    public static float SHBasisL0()
    public static float SHBasisL1_1(Vector3 direction)
    public static float SHBasisL10(Vector3 direction)
    public static float SHBasisL11(Vector3 direction)
    public static float SHBasisL2_2(Vector3 direction)
    public static float SHBasisL2_1(Vector3 direction)
    public static float SHBasisL20(Vector3 direction)
    public static float SHBasisL21(Vector3 direction)
    public static float SHBasisL22(Vector3 direction)
    public static float SHBasis(int index, Vector3 direction)
    
    public static void UnityConventionToRawCoefficientsInPlace(ref SphericalHarmonicsL2 sh)
    public static void RawCoefficientsToUnityConventionInPlace(ref SphericalHarmonicsL2 sh)
    public static SphericalHarmonicsL2 UnityConventionToRawCoefficients(in SphericalHarmonicsL2 sh)
    public static SphericalHarmonicsL2 RawCoefficientsToUnityConvention(in SphericalHarmonicsL2 sh)
    
    public static SphericalHarmonicsL2 Lerp(in SphericalHarmonicsL2 a, in SphericalHarmonicsL2 b, float t)
    
    public static Vector4[] GetShaderCoefficients(in SphericalHarmonicsL2 sh)
}

struct RawSphericalHarmonicsL2 {
    public readonly SphericalHarmonicsL2 AsUnityConvention()
    public readonly SphericalHarmonicsL2 AsRaw()
    public Color Evaluate(Vector3 direction)
    public void Clear()

    public static RawSphericalHarmonicsL2 Lerp(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b, float t)

    public static RawSphericalHarmonicsL2 Convolve(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b);
    public void Convolve(in RawSphericalHarmonicsL2 other);
    public static RawSphericalHarmonicsL2 Product(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b);
    public void Product(in RawSphericalHarmonicsL2 other);

    public static void ConvolveRadianceToIrradianceInPlace(ref RawSphericalHarmonicsL2 sh)
    public static RawSphericalHarmonicsL2 ConvolveRadianceToIrradiance(in RawSphericalHarmonicsL2 sh)
    public static void DeConvolveIrradianceToRadianceInPlace(ref RawSphericalHarmonicsL2 sh)
    public static RawSphericalHarmonicsL2 DeConvolveIrradianceToRadiance(in RawSphericalHarmonicsL2 sh)

    public static RawSphericalHarmonicsL2 ProjectIntoSHMonteCarlo(Func<Vector3, Color> sphericalFunction, int sampleCount)
    public static RawSphericalHarmonicsL2 ProjectIntoSHMonteCarlo(Func<Vector3, Color> sphericalFunction, Func<int, Vector3> rngFunction, int sampleCount)
    public static RawSphericalHarmonicsL2 ProjectIntoSHRiemann(Func<Vector3, Color> sphericalFunction, int samplesPhi, int samplesTheta)

    public static RawSphericalHarmonicsL2 ProjectCubemapIntoSHMonteCarlo(Cubemap cubemap, int sampleCount, bool convolveToIrradiance)
    public static RawSphericalHarmonicsL2 ProjectCubemapIntoSHMonteCarlo(Cubemap cubemap, Func<int, Vector3> rngFunction, int sampleCount, bool convolveToIrradiance)
    public static RawSphericalHarmonicsL2 ProjectCubemapIntoSHRiemann(Cubemap cubemap, int samplesPhi, int samplesTheta, bool convolveToIrradiance)
}

struct SHMatrix {
    public static SHMatrix RotateZ(float angle)
    public static readonly SHMatrix RotationXPositive90;
    public static readonly SHMatrix RotationXNegative90;

    public static SHMatrix RotateZYZ(Vector3 alphaBetaGammaRadians)
    public static SHMatrix Rotate(Vector3 xyzDegrees)
}
```

# Examples
Rotating a light probe:
```cs
SphericalHarmonicsL2 sh = LightmapSettings.lightProbes.bakedProbes[0];
SHMatrix rotMat = SHMatrix.Rotate(transform.eulerAngles);
sh = rotMat * sh;
LightmapSettings.lightProbes.bakedProbes = new[] { sh };
```

Projecting a function into SH:
```cs
SphericalHarmonicsL2 sh1 = RawSphericalHarmonicsL2.ProjectIntoSHMonteCarlo(
    (direction) => {
        Vector3 dirScaled = direction * 0.5f + new Vector3(0.5f, 0.5f, 0.5f);
        return new Color(dirScaled.x, dirScaled.y, dirScaled.z);
    },
    1000);

SphericalHarmonicsL2 sh2 = RawSphericalHarmonicsL2.ProjectIntoSHRiemann(
    (direction) => {
        Vector3 dirScaled = direction * 0.5f + new Vector3(0.5f, 0.5f, 0.5f);
        return new Color(dirScaled.x, dirScaled.y, dirScaled.z);
    },
    100, 100);
```
