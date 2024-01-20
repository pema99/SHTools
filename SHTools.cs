using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Unity.Collections.LowLevel.Unsafe;
using UnityEngine;
using UnityEngine.Rendering;

namespace SHTools
{
    public static class SHUtility
    {
        #region SH Basis functions
        /// <summary>
        /// Constant part of the SH basis function Y(0,0) = 1/2 * sqrt(1/π).
        /// </summary>
        public const float SH_L0_NORMALIZATION = 0.2820947917738781434740397257803862929220253146644994284220428608f;
        /// <summary>
        /// Constant part of the SH basis function Y(1,x) = 1/2 * sqrt(3/π).
        /// </summary>
        public const float SH_L1_NORMALIZATION = 0.4886025119029199215863846228383470045758856081942277021382431574f;
        /// <summary>
        /// Constant part of the SH basis function Y(2,-2) = sqrt(15/π)/2.
        /// </summary>
        public const float SH_L2_2_NORMALIZATION = 1.0925484305920790705433857058026884026904329595042589753478516999f;
        /// <summary>
        /// Constant part of the SH basis function Y(2,-1) = sqrt(15/π)/2.
        /// </summary>
        public const float SH_L2_1_NORMALIZATION = SH_L2_2_NORMALIZATION;
        /// <summary>
        /// Constant part of the SH basis function Y(2,0) = sqrt(5/π)/4.
        /// </summary>
        public const float SH_L20_NORMALIZATION = 0.3153915652525200060308936902957104933242475070484115878434078878f;
        /// <summary>
        /// Constant part of the SH basis function Y(2,1) = sqrt(15/π)/2.
        /// </summary>
        public const float SH_L21_NORMALIZATION = SH_L2_2_NORMALIZATION;
        /// <summary>
        /// Constant part of the SH basis function Y(2,2) = sqrt(15/π)/4.
        /// </summary>
        public const float SH_L22_NORMALIZATION = 0.5462742152960395352716928529013442013452164797521294876739258499f;

        /// <summary>
        /// SH Basis function Y(0,0). This is the constant term.
        /// </summary>
        /// <returns>Y(0,0) evaluted.</returns>
        public static float SHBasisL0()
        {
            return SH_L0_NORMALIZATION;
        }

        /// <summary>
        /// SH Basis function Y(1,-1).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(1,-1) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL1_1(Vector3 direction)
        {
            return SH_L1_NORMALIZATION * direction.x;
        }

        /// <summary>
        /// SH Basis function Y(1,0).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(1,0) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL10(Vector3 direction)
        {
            return SH_L1_NORMALIZATION * direction.y;
        }

        /// <summary>
        /// SH Basis function Y(1,1).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(1,1) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL11(Vector3 direction)
        {
            return SH_L1_NORMALIZATION * direction.z;
        }

        /// <summary>
        /// SH Basis function Y(2,-2).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(2,-2) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL2_2(Vector3 direction)
        {
            return SH_L2_2_NORMALIZATION * direction.x * direction.y;
        }

        /// <summary>
        /// SH Basis function Y(2,-1).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(2,-1) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL2_1(Vector3 direction)
        {
            return SH_L2_1_NORMALIZATION * direction.y * direction.z;
        }

        /// <summary>
        /// SH Basis function Y(2,0).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(2,0) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL20(Vector3 direction)
        {
            return SH_L20_NORMALIZATION * (3.0f * direction.z * direction.z - 1.0f);
        }

        /// <summary>
        /// SH Basis function Y(2,1).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(2,1) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL21(Vector3 direction)
        {
            return SH_L21_NORMALIZATION * direction.x * direction.z;
        }

        /// <summary>
        /// SH Basis function Y(2,2).
        /// </summary>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Y(2,2) evaluated at the given <paramref name="direction"/>.</returns>
        public static float SHBasisL22(Vector3 direction)
        {
            return SH_L22_NORMALIZATION * (direction.x * direction.x - direction.y * direction.y);
        } 

        /// <summary>
        /// SH Basis function at the given <paramref name="index"/>.
        /// </summary>
        /// <param name="index">Index of the SH basis function.</param>
        /// <param name="direction">Direction to evaluate the SH basis function at.</param>	 
        /// <returns>Basis function at the given <paramref name="index"/> evaluated at the given <paramref name="direction"/>.</returns>
        /// <exception cref="ArgumentException">Thrown if the given <paramref name="index"/> isn't in range [0; 8].</exception>
        public static float SHBasis(int index, Vector3 direction)
        {
            switch (index)
            {
                case 0: return SHBasisL0();
                case 1: return SHBasisL1_1(direction);
                case 2: return SHBasisL10(direction);
                case 3: return SHBasisL11(direction);
                case 4: return SHBasisL2_2(direction);
                case 5: return SHBasisL2_1(direction);
                case 6: return SHBasisL20(direction);
                case 7: return SHBasisL21(direction);
                case 8: return SHBasisL22(direction);
                default: throw new ArgumentException("Invalid SH basis index");
            }
        }
        #endregion

        #region General utilities
        /// <summary>
        /// Unity doesn't store raw SH coefficients in SphericalHarmonicsL2. Instead, it stores the
        /// coefficients multiplied by the constant part of each basis function, and divided by PI.
        /// This function converts from Unity's convention back to the raw coefficients.
        /// </summary>
        /// <param name="sh">The SH to convert.</param>
        public static void UnityConventionToRawCoefficientsInPlace(ref SphericalHarmonicsL2 sh)
        {
            for (int i = 0; i < 3; i++)
            {
                sh[i, 0] = (sh[i, 0] * Mathf.PI) / SHUtility.SH_L0_NORMALIZATION;

                sh[i, 1] = (sh[i, 1] * Mathf.PI) / SHUtility.SH_L1_NORMALIZATION;
                sh[i, 2] = (sh[i, 2] * Mathf.PI) / SHUtility.SH_L1_NORMALIZATION;
                sh[i, 3] = (sh[i, 3] * Mathf.PI) / SHUtility.SH_L1_NORMALIZATION;

                sh[i, 4] = (sh[i, 4] * Mathf.PI) / SHUtility.SH_L2_2_NORMALIZATION;
                sh[i, 5] = (sh[i, 5] * Mathf.PI) / SHUtility.SH_L2_1_NORMALIZATION;
                sh[i, 6] = (sh[i, 6] * Mathf.PI) / SHUtility.SH_L20_NORMALIZATION;
                sh[i, 7] = (sh[i, 7] * Mathf.PI) / SHUtility.SH_L21_NORMALIZATION;
                sh[i, 8] = (sh[i, 8] * Mathf.PI) / SHUtility.SH_L22_NORMALIZATION;
            }
        }

        /// <summary>
        /// Unity doesn't store raw SH coefficients in SphericalHarmonicsL2. Instead, it stores the
        /// coefficients multiplied by the constant part of each basis function, and divided by PI.
        /// This function converts from raw coefficients to Unity's convention.
        /// </summary>
        /// <param name="sh">The SH to convert.</param>
        public static void RawCoefficientsToUnityConventionInPlace(ref SphericalHarmonicsL2 sh)
        {
            for (int i = 0; i < 3; i++)
            {
                sh[i, 0] = (sh[i, 0] * SHUtility.SH_L0_NORMALIZATION) / Mathf.PI;

                sh[i, 1] = (sh[i, 1] * SHUtility.SH_L1_NORMALIZATION) / Mathf.PI;
                sh[i, 2] = (sh[i, 2] * SHUtility.SH_L1_NORMALIZATION) / Mathf.PI;
                sh[i, 3] = (sh[i, 3] * SHUtility.SH_L1_NORMALIZATION) / Mathf.PI;

                sh[i, 4] = (sh[i, 4] * SHUtility.SH_L2_2_NORMALIZATION) / Mathf.PI;
                sh[i, 5] = (sh[i, 5] * SHUtility.SH_L2_1_NORMALIZATION) / Mathf.PI;
                sh[i, 6] = (sh[i, 6] * SHUtility.SH_L20_NORMALIZATION) / Mathf.PI;
                sh[i, 7] = (sh[i, 7] * SHUtility.SH_L21_NORMALIZATION) / Mathf.PI;
                sh[i, 8] = (sh[i, 8] * SHUtility.SH_L22_NORMALIZATION) / Mathf.PI;
            }
        }

        /// <summary>
        /// Unity doesn't store raw SH coefficients in SphericalHarmonicsL2. Instead, it stores the
        /// coefficients multiplied by the constant part of each basis function, and divided by PI.
        /// This function converts from Unity's convention back to the raw coefficients.  
        /// </summary>
        /// <param name="sh">The SH to convert.</param>
        public static SphericalHarmonicsL2 UnityConventionToRawCoefficients(in SphericalHarmonicsL2 sh)
        {
            SphericalHarmonicsL2 result = sh;
            UnityConventionToRawCoefficientsInPlace(ref result);
            return result;
        }

        /// <summary>
        /// Unity doesn't store raw SH coefficients in SphericalHarmonicsL2. Instead, it stores the
        /// coefficients multiplied by the constant part of each basis function, and divided by PI.
        /// This function converts from raw coefficients to Unity's convention.  
        /// </summary>
        /// <param name="sh">The SH to convert.</param>
        public static SphericalHarmonicsL2 RawCoefficientsToUnityConvention(in SphericalHarmonicsL2 sh)
        {
            SphericalHarmonicsL2 result = sh;
            RawCoefficientsToUnityConventionInPlace(ref result);
            return result;
        }

        /// <summary>
        /// Linearly interpolate between two sets of SH coefficients.
        /// </summary>
        /// <param name="a">First SH.</param>
        /// <param name="b">Second SH.</param>
        /// <param name="t">Interpolation parameter.</param>
        /// <returns>Interpolated set of SH coefficients.</returns>
        public static SphericalHarmonicsL2 Lerp(in SphericalHarmonicsL2 a, in SphericalHarmonicsL2 b, float t)
        {
            SphericalHarmonicsL2 result = new();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    result[i, j] = Mathf.Lerp(a[i, j], b[i, j], t);
                }
            }
            return result;
        }

        /// <summary>
        /// Get's the shader coefficients (unity_SHAr...unity_SHC) for this SphericalHarmonicsL2.
        /// </summary>
        /// <param name="sh">SH to get shader coefficients for.</param>
        /// <returns>Shader coefficients in order unity_SHAr, unity_SHAg, unity_SHAb, unity_SHBr, unity_SHBg, unity_SHBb, unity_SHC.</returns> 
        public static Vector4[] GetShaderCoefficients(in SphericalHarmonicsL2 sh)
        {
            Vector4[] result = new Vector4[7];
            for (int i = 0; i < 3; i++)
            {
                result[i] = new Vector4(
                    sh[i, 3],
                    sh[i, 1],
                    sh[i, 2],
                    sh[i, 0] - sh[i, 6]
                );

                result[i + 3] = new Vector4(
                    sh[i, 4],
                    sh[i, 5],
                    sh[i, 6] * 3.0f,
                    sh[i, 7]
                );
            }

            result[6] = new Vector4(
                sh[0, 8],
                sh[1, 8],
                sh[2, 8],
                1.0f
            );
            return result;
        }

        /// <summary>
        /// Project a spherical function into SH using Monte Carlo integration.
        /// </summary>
        /// <param name="sphericalFunction">Function to project, going from a direction to a color.</param>
        /// <param name="sampleCount">Number of samples to use.</param>
        /// <returns>The function projected into a RawSphericalHarmonicsL2.</returns>
        public static RawSphericalHarmonicsL2 ProjectIntoSHMonteCarlo(Func<Vector3, Vector3> sphericalFunction, int sampleCount)
        {
            return ProjectIntoSHMonteCarlo(sphericalFunction, i => UnityEngine.Random.onUnitSphere, sampleCount);
        }

        /// <summary>
        /// Project a spherical function into SH using Monte Carlo integration.
        /// </summary>
        /// <param name="sphericalFunction">Function to project, going from a direction to a color.</param>
        /// <param name="rngFunction">Function to generate random direction vectors, given the sample index.</param>
        /// <param name="sampleCount">Number of samples to use.</param>
        /// <returns>The function projected into a RawSphericalHarmonicsL2.</returns>
        public static RawSphericalHarmonicsL2 ProjectIntoSHMonteCarlo(Func<Vector3, Vector3> sphericalFunction, Func<int, Vector3> rngFunction, int sampleCount)
        {
            RawSphericalHarmonicsL2 result = new();

            for (int i = 0; i < sampleCount; i++)
            {
                Vector3 direction = rngFunction(i);
                Vector3 eval = sphericalFunction(direction);

                for (int n = 0; n < 9; n++)
                {
                    Vector3 proj = eval * SHBasis(n, direction);
                    result[0, n] += proj.x;
                    result[1, n] += proj.y;
                    result[2, n] += proj.z;
                }
            }
            
            // Monte carlo normalization
            float reciprocalSampleCount = 1.0f / sampleCount;
            float reciprocalUniformSphereDensity = 4.0f * Mathf.PI;
            return result * reciprocalSampleCount * reciprocalUniformSphereDensity;
        }
        #endregion
    }

    /// <summary>
    /// Unity doesn't store raw SH coefficients in SphericalHarmonicsL2. Instead, it stores the
    /// coefficients multiplied by the constant part of each basis function, and divided by PI.
    /// Crucially, SH coefficients in Unity's convention cannot directly be rotated. 
    /// This type represents a set of raw SH coefficients, which can be rotated, and may freely
    /// be converted to and from Unity's convention. 
    /// </summary>
    public struct RawSphericalHarmonicsL2
    {
        private SphericalHarmonicsL2 sh;

        /// <summary>
        /// Construct a RawSphericalHarmonicsL2 from a SphericalHarmonicsL2 in Unity's convention.
        /// </summary>
        /// <param name="unityConventionSH">SH using Unity's convention.</param>
        public RawSphericalHarmonicsL2(in SphericalHarmonicsL2 unityConventionSH)
        {
            sh = SHUtility.UnityConventionToRawCoefficients(unityConventionSH);
        }

        public static implicit operator SphericalHarmonicsL2(RawSphericalHarmonicsL2 rawSH)
        {
            return SHUtility.RawCoefficientsToUnityConvention(rawSH.sh);
        }

        public static implicit operator RawSphericalHarmonicsL2(SphericalHarmonicsL2 unityConventionSH)
        {
            return new RawSphericalHarmonicsL2(unityConventionSH);
        }

        /// <summary>
        /// Get the SH coefficients in Unity's convention.
        /// </summary>
        /// <returns>SH coefficients in Unity's convention.</returns>
        public readonly SphericalHarmonicsL2 AsUnityConvention()
        {
            return SHUtility.RawCoefficientsToUnityConvention(sh);
        }

        /// <summary>
        /// Get the raw underlying SH coefficients.
        /// </summary>
        /// <returns>Raw SH coefficients.</returns>
        public readonly SphericalHarmonicsL2 AsRaw()
        {
            return sh;
        }

        /// <summary>
        /// Get or set the coefficient at the given <paramref name="channel"/> and <paramref name="coefficient"/>.
        /// </summary>
        public float this[int channel, int coefficient]
        {
            get => sh[channel, coefficient];
            set => sh[channel, coefficient] = value;
        }

        /// <summary>
        /// Linearly interpolate between two sets of SH coefficients.
        /// </summary>
        /// <param name="a">First SH.</param>
        /// <param name="b">Second SH.</param>
        /// <param name="t">Interpolation parameter.</param>
        /// <returns>Interpolated set of SH coefficients.</returns>
        public static RawSphericalHarmonicsL2 Lerp(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b, float t)
        {
            return SHUtility.Lerp(a.AsRaw(), b.AsRaw(), t);
        }

        // Radiance convolution constants
        private const float AHat0 = Mathf.PI;
        private const float AHat1 = 2.09439510239319549230842892218633525613f; // 2*PI/3
        private const float AHat2 = 0.78539816339744830961566084581987572104f; // PI/4

        /// <summary>
        /// Convolves radiance stored in SH to irradiance, in place.
        /// Based on https://cseweb.ucsd.edu/~ravir/papers/envmap/envmap.pdf
        /// </summary>
        /// <param name="sh">SH to convolve.</param>
        public static void ConvolveRadianceToIrradianceInPlace(ref RawSphericalHarmonicsL2 sh)
        {
            for (int i = 0; i < 3; i++)
            {
                sh[i, 0] *= AHat0;
                sh[i, 1] *= AHat1;
                sh[i, 2] *= AHat1;
                sh[i, 3] *= AHat1;
                sh[i, 4] *= AHat2;
                sh[i, 5] *= AHat2;
                sh[i, 6] *= AHat2;
                sh[i, 7] *= AHat2;
                sh[i, 8] *= AHat2;
            }
        }

        /// <summary>
        /// Convolves radiance stored in SH to irradiance.
        /// Based on https://cseweb.ucsd.edu/~ravir/papers/envmap/envmap.pdf
        /// </summary>
        /// <param name="sh">SH to convolve.</param>
        public static RawSphericalHarmonicsL2 ConvolveRadianceToIrradiance(in RawSphericalHarmonicsL2 sh)
        {
            RawSphericalHarmonicsL2 result = sh;
            ConvolveRadianceToIrradianceInPlace(ref result);
            return result;
        }

        /// <summary>
        /// Converts irradiance stored in SH back to radiance, in place.
        /// Based on https://cseweb.ucsd.edu/~ravir/papers/envmap/envmap.pdf
        /// </summary>
        /// <param name="sh">SH to deconvolve.</param>
        public static void DeConvolveIrradianceToRadianceInPlace(ref RawSphericalHarmonicsL2 sh)
        {
            for (int i = 0; i < 3; i++)
            {
                sh[i, 0] /= AHat0;
                sh[i, 1] /= AHat1;
                sh[i, 2] /= AHat1;
                sh[i, 3] /= AHat1;
                sh[i, 4] /= AHat2;
                sh[i, 5] /= AHat2;
                sh[i, 6] /= AHat2;
                sh[i, 7] /= AHat2;
                sh[i, 8] /= AHat2;
            }
        }

        /// <summary>
        /// Converts irradiance stored in SH back to radiance.
        /// Based on https://cseweb.ucsd.edu/~ravir/papers/envmap/envmap.pdf
        /// </summary>
        /// <param name="sh">SH to deconvolve.</param>
        public static RawSphericalHarmonicsL2 DeConvolveIrradianceToRadiance(in RawSphericalHarmonicsL2 sh)
        {
            RawSphericalHarmonicsL2 result = sh;
            DeConvolveIrradianceToRadianceInPlace(ref result);
            return result;
        }

        /// <summary>
        /// Evaluate the SH at the given <paramref name="direction"/>.
        /// </summary>
        /// <param name="direction">Direction to evaluate at.</param>
        /// <returns>Color result of evaluation.</returns>
        public Color Evaluate(Vector3 direction)
        {
            Color result = new(0,0,0);
            for (int channel = 0; channel < 3; channel++)
            {
                for (int index = 0; index < 9; index++)
                {
                    float coefficient = sh[channel, index];
                    float basis = SHUtility.SHBasis(index, direction);
                    result[channel] += coefficient * basis;
                }
            }
            return result;
        }

        /// <summary>
        /// Clears the SH coefficients to zero.
        /// </summary>
        public void Clear()
        {
            sh.Clear();
        }

        public static RawSphericalHarmonicsL2 operator +(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b)
        {
            RawSphericalHarmonicsL2 result = new();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    result[i, j] = a[i, j] + b[i, j];
                }
            }
            return result;
        }

        public static RawSphericalHarmonicsL2 operator *(in RawSphericalHarmonicsL2 a, float b)
        {
            RawSphericalHarmonicsL2 result = new();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    result[i, j] = a[i, j] * b;
                }
            }
            return result;
        }

        public static bool operator ==(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b)
        {
            return a.sh == b.sh;
        }

        public static bool operator !=(in RawSphericalHarmonicsL2 a, in RawSphericalHarmonicsL2 b)
        {
            return a.sh != b.sh;
        }

        public override readonly bool Equals(object obj)
        {
            return obj is RawSphericalHarmonicsL2 sh && this == sh;
        }

        public override int GetHashCode()
        {
            return sh.GetHashCode();
        }
    }

    /// <summary>
    /// SH basis levels, supports up to L2.
    /// </summary>
    public enum SHLevel
    {
        L0,
        L1,
        L2
    }

    /// <summary>
    /// Represents a transformation on L2 spherical harmonics in block matrix form.
    /// Consists of 2 blocks: The L1 and L2 block. L0 is never transformed. 
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct SHMatrix
    {
        // L1 block
        public float L1M00, L1M01, L1M02;
        public float L1M10, L1M11, L1M12;
        public float L1M20, L1M21, L1M22;

        // L2 block
        public float L2M00, L2M01, L2M02, L2M03, L2M04;
        public float L2M10, L2M11, L2M12, L2M13, L2M14;
        public float L2M20, L2M21, L2M22, L2M23, L2M24;
        public float L2M30, L2M31, L2M32, L2M33, L2M34;
        public float L2M40, L2M41, L2M42, L2M43, L2M44;

        private float GetL1Coefficient(int row, int col)
        {
            switch (row * 3 + col)
            {
                case 0: return L1M00;
                case 1: return L1M01;
                case 2: return L1M02;
                case 3: return L1M10;
                case 4: return L1M11;
                case 5: return L1M12;
                case 6: return L1M20;
                case 7: return L1M21;
                case 8: return L1M22;
                default: throw new ArgumentException("Invalid SH level");
            }
        }

        private void SetL1Coefficient(int row, int col, float value)
        {
            switch (row * 3 + col)
            {
                case 0: L1M00 = value; break;
                case 1: L1M01 = value; break;
                case 2: L1M02 = value; break;
                case 3: L1M10 = value; break;
                case 4: L1M11 = value; break;
                case 5: L1M12 = value; break;
                case 6: L1M20 = value; break;
                case 7: L1M21 = value; break;
                case 8: L1M22 = value; break;
                default: throw new ArgumentException("Invalid SH level");
            }
        }

        private float GetL2Coefficient(int row, int col)
        {
            switch (row * 5 + col)
            {
                case 0: return L2M00;
                case 1: return L2M01;
                case 2: return L2M02;
                case 3: return L2M03;
                case 4: return L2M04;
                case 5: return L2M10;
                case 6: return L2M11;
                case 7: return L2M12;
                case 8: return L2M13;
                case 9: return L2M14;
                case 10: return L2M20;
                case 11: return L2M21;
                case 12: return L2M22;
                case 13: return L2M23;
                case 14: return L2M24;
                case 15: return L2M30;
                case 16: return L2M31;
                case 17: return L2M32;
                case 18: return L2M33;
                case 19: return L2M34;
                case 20: return L2M40;
                case 21: return L2M41;
                case 22: return L2M42;
                case 23: return L2M43;
                case 24: return L2M44;
                default: throw new ArgumentException("Invalid SH level");
            }
        }

        private void SetL2Coefficient(int row, int col, float value)
        {
            switch (row * 5 + col)
            {
                case 0: L2M00 = value; break;
                case 1: L2M01 = value; break;
                case 2: L2M02 = value; break;
                case 3: L2M03 = value; break;
                case 4: L2M04 = value; break;
                case 5: L2M10 = value; break;
                case 6: L2M11 = value; break;
                case 7: L2M12 = value; break;
                case 8: L2M13 = value; break;
                case 9: L2M14 = value; break;
                case 10: L2M20 = value; break;
                case 11: L2M21 = value; break;
                case 12: L2M22 = value; break;
                case 13: L2M23 = value; break;
                case 14: L2M24 = value; break;
                case 15: L2M30 = value; break;
                case 16: L2M31 = value; break;
                case 17: L2M32 = value; break;
                case 18: L2M33 = value; break;
                case 19: L2M34 = value; break;
                case 20: L2M40 = value; break;
                case 21: L2M41 = value; break;
                case 22: L2M42 = value; break;
                case 23: L2M43 = value; break;
                case 24: L2M44 = value; break;
                default: throw new ArgumentException("Invalid SH level");
            }
        }

        /// <summary>
        /// Get or set the coefficient at the given <paramref name="level"/>, <paramref name="row"/> and <paramref name="col"/>.
        /// </summary>
        public float this[SHLevel level, int row, int col]
        {
            get
            {
                switch (level)
                {
                    case SHLevel.L0: return 1.0f;
                    case SHLevel.L1: return GetL1Coefficient(row, col);
                    case SHLevel.L2: return GetL2Coefficient(row, col);
                    default: throw new ArgumentException("Invalid SH level");
                }
            }
            set
            {
                switch (level)
                {
                    case SHLevel.L0: throw new ArgumentException("Cannot set L0");
                    case SHLevel.L1: SetL1Coefficient(row, col, value); break;
                    case SHLevel.L2: SetL2Coefficient(row, col, value); break;
                    default: throw new ArgumentException("Invalid SH level");
                }
            }
        }

        private static void MultiplyBlock(SHLevel level, in SHMatrix a, in SHMatrix b, ref SHMatrix result)
        {
            if (level == SHLevel.L0) return;
            int size = level == SHLevel.L1 ? 3 : 5;

            for (int row = 0; row < size; row++)
            {
                for (int col = 0; col < size; col++) 
                {
                    float sum = 0.0f;
                    for (int i = 0; i < size; i++)
                    {
                        sum += a[level, row, i] * b[level, i, col];
                    }
                    result[level, row, col] = sum;
                }
            }
        }

        private static void MultiplyBlock(SHLevel level, in SHMatrix a, in RawSphericalHarmonicsL2 b, ref RawSphericalHarmonicsL2 result)
        {
            if (level == SHLevel.L0) return;
            int size = level == SHLevel.L1 ? 3 : 5;
            int offset = level == SHLevel.L1 ? 1 : 4;

            for (int channel = 0; channel < 3; channel++)
            {
                for (int row = 0; row < size; row++)
                {
                    float sum = 0.0f;
                    for (int col = 0; col < size; col++)
                    {
                        sum += a[level, row, col] * b[channel, col + offset];
                    }
                    result[channel, row + offset] = sum;
                }
            }
        }

        public static SHMatrix operator *(in SHMatrix a, in SHMatrix b)
        {
            SHMatrix result = new();
            MultiplyBlock(SHLevel.L1, in a, in b, ref result);
            MultiplyBlock(SHLevel.L2, in a, in b, ref result);
            return result;
        }

        public static RawSphericalHarmonicsL2 operator *(in SHMatrix a, in RawSphericalHarmonicsL2 b)
        {
            RawSphericalHarmonicsL2 result = new();
            result[0, 0] = b[0, 0];
            result[1, 0] = b[1, 0];
            result[2, 0] = b[2, 0];
            MultiplyBlock(SHLevel.L1, in a, in b, ref result);
            MultiplyBlock(SHLevel.L2, in a, in b, ref result);
            return result;
        }

        /// <summary>
        /// Creates a rotation matrix around the Z axis.
        /// </summary>
        /// <param name="angle">Angle to rotate in radians.</param>
        /// <returns>Rotation matrix representing rotation around the Z axis.</returns>
        public static SHMatrix RotateZ(float angle)
        {
            SHMatrix matrix = new();

            float cosA = Mathf.Cos(angle);
            float sinA = Mathf.Sin(angle);
            float cos2A = Mathf.Cos(2.0f * angle);
            float sin2A = Mathf.Sin(2.0f * angle);

            matrix.L1M00 = cosA;  matrix.L1M01 = 0; matrix.L1M02 = sinA;
            matrix.L1M10 = 0;     matrix.L1M11 = 1; matrix.L1M12 = 0;
            matrix.L1M20 = -sinA; matrix.L1M21 = 0; matrix.L1M22 = cosA;

            matrix.L2M00 = cos2A;  matrix.L2M01 = 0;     matrix.L2M02 = 0; matrix.L2M03 = 0;    matrix.L2M04 = sin2A;
            matrix.L2M10 = 0;      matrix.L2M11 = cosA;  matrix.L2M12 = 0; matrix.L2M13 = sinA; matrix.L2M14 = 0;
            matrix.L2M20 = 0;      matrix.L2M21 = 0;     matrix.L2M22 = 1; matrix.L2M23 = 0;    matrix.L2M24 = 0;
            matrix.L2M30 = 0;      matrix.L2M31 = -sinA; matrix.L2M32 = 0; matrix.L2M33 = cosA; matrix.L2M34 = 0;
            matrix.L2M40 = -sin2A; matrix.L2M41 = 0;     matrix.L2M42 = 0; matrix.L2M43 = 0;    matrix.L2M44 = cos2A;

            return matrix;
        }

        /// <summary>
        /// The rotation matrix for a 90 degree rotation around the X axis.
        /// </summary>
        public static readonly SHMatrix RotationXPositive90 = new()
        {
            L1M00 = 0, L1M01 = -1, L1M02 = 0,
            L1M10 = 1, L1M11 = 0,  L1M12 = 0,
            L1M20 = 0, L1M21 = 0,  L1M22 = 1,

            L2M00 = 0, L2M01 = 0,  L2M02 = 0,     L2M03 = -1, L2M04 = 0,
            L2M10 = 0, L2M11 = -1, L2M12 = 0,     L2M13 = 0,  L2M14 = 0,
            L2M20 = 0, L2M21 = 0,  L2M22 = -0.5f, L2M23 = 0,  L2M24 = -Mathf.Sqrt(3.0f) / 2.0f,
            L2M30 = 1, L2M31 = 0,  L2M32 = 0,     L2M33 = 0,  L2M34 = 0,
            L2M40 = 0, L2M41 = 0,  L2M42 = -Mathf.Sqrt(3.0f) / 2.0f, L2M43 = 0, L2M44 = 0.5f
        };

        /// <summary>
        /// The rotation matrix for a -90 degree rotation around the X axis.
        /// </summary>
        public static readonly SHMatrix RotationXNegative90 = new()
        {
            L1M00 = 0,  L1M01 = 1, L1M02 = 0,
            L1M10 = -1, L1M11 = 0, L1M12 = 0,
            L1M20 = 0,  L1M21 = 0, L1M22 = 1,

            L2M00 = 0, L2M01 = 0,  L2M02 = 0,     L2M03 = 1, L2M04 = 0,
            L2M10 = 0, L2M11 = -1, L2M12 = 0,     L2M13 = 0, L2M14 = 0,
            L2M20 = 0, L2M21 = 0,  L2M22 = -0.5f, L2M23 = 0, L2M24 = -Mathf.Sqrt(3.0f) / 2.0f,
            L2M30 = -1, L2M31 = 0, L2M32 = 0,    L2M33 = 0, L2M34 = 0,
            L2M40 = 0, L2M41 = 0,  L2M42 = -Mathf.Sqrt(3.0f) / 2.0f, L2M43 = 0, L2M44 = 0.5f
        };

        /// <summary>
        /// Creates a rotation matrix from the given euler angles in radians using ZYZ convention.
        /// </summary>
        /// <param name="alphaBetaGamma">Euler angles in radians using ZYZ convention.</param>
        /// <returns>Rotation matrix.</returns>
        public static SHMatrix RotateZYZ(Vector3 alphaBetaGammaRadians)
        {
            return
                RotateZ(alphaBetaGammaRadians.z) *
                RotationXNegative90 *
                RotateZ(alphaBetaGammaRadians.y) *
                RotationXPositive90 *
                RotateZ(alphaBetaGammaRadians.x);
        }

        private static Vector3 EulerAnglesToZYZ(Vector3 xyzDegrees)
        {
            var eulerMat = Matrix4x4.Rotate(Quaternion.Euler(xyzDegrees.x, xyzDegrees.y, xyzDegrees.z));

            float cosBeta = eulerMat[2, 2];
            float sinBeta = Mathf.Sqrt(1.0f - eulerMat[2, 2] * eulerMat[2, 2]);
            float beta = -Mathf.Atan2(sinBeta, cosBeta);
            
            if (sinBeta < 0.001f)
            {
                float cosAlpha = eulerMat[1, 1];
                float sinAlpha = -eulerMat[1, 0];
                float alpha = -Mathf.Atan2(sinAlpha, cosAlpha);

                float cosGamma = 1.0f;
                float sinGamma = 0.0f;
                float gamma = -Mathf.Atan2(sinGamma, cosGamma);

                return new Vector3(alpha, beta, gamma);
            }
            else
            {
                float cosAlpha = eulerMat[2, 0] / sinBeta;
                float sinAlpha = eulerMat[2, 1] / sinBeta;
                float alpha = -Mathf.Atan2(sinAlpha, cosAlpha);

                float cosGamma = -eulerMat[0, 2] / sinBeta;
                float sinGamma = eulerMat[1, 2] / sinBeta;
                float gamma = -Mathf.Atan2(sinGamma, cosGamma);

                return new Vector3(alpha, beta, gamma);
            }
        }

        /// <summary>
        /// Creates a rotation matrix from the given euler angles in degrees using XYZ convention.
        /// </summary>
        /// <param name="alphaBetaGamma">Euler angles in degrees using XYZ convention.</param>
        /// <returns>Rotation matrix.</returns>
        public static SHMatrix Rotate(Vector3 xyzDegrees)
        {
            Vector3 alphaBetaGamma = EulerAnglesToZYZ(xyzDegrees);
            return RotateZYZ(alphaBetaGamma);
        }
    }
}