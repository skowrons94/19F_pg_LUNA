import numba
import numpy as np
from scipy import constants
from numba.experimental import jitclass

from pynucastro.rates import TableIndex, TableInterpolator, TabularRate, Tfactors
from pynucastro.screening import PlasmaState, ScreenFactors

jp = 0
jhe4 = 1
jo16 = 2
jo17 = 3
jo18 = 4
jf17 = 5
jf18 = 6
jf19 = 7
jne20 = 8
jne21 = 9
jne22 = 10
jna21 = 11
jna22 = 12
jna23 = 13
nnuc = 14

A = np.zeros((nnuc), dtype=np.int32)

A[jp] = 1
A[jhe4] = 4
A[jo16] = 16
A[jo17] = 17
A[jo18] = 18
A[jf17] = 17
A[jf18] = 18
A[jf19] = 19
A[jne20] = 20
A[jne21] = 21
A[jne22] = 22
A[jna21] = 21
A[jna22] = 22
A[jna23] = 23

Z = np.zeros((nnuc), dtype=np.int32)

Z[jp] = 1
Z[jhe4] = 2
Z[jo16] = 8
Z[jo17] = 8
Z[jo18] = 8
Z[jf17] = 9
Z[jf18] = 9
Z[jf19] = 9
Z[jne20] = 10
Z[jne21] = 10
Z[jne22] = 10
Z[jna21] = 11
Z[jna22] = 11
Z[jna23] = 11

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jp] = 0.0015040963030260536
mass[jhe4] = 0.0059735574925878256
mass[jo16] = 0.023871099858982767
mass[jo17] = 0.02536981167252093
mass[jo18] = 0.02686227133140636
mass[jf17] = 0.025374234423440733
mass[jf18] = 0.026864924401329426
mass[jf19] = 0.028353560468882166
mass[jne20] = 0.02983707929641827
mass[jne21] = 0.03133159647374143
mass[jne22] = 0.0328203408644564
mass[jna21] = 0.0313372792660881
mass[jna22] = 0.03282489638134515
mass[jna23] = 0.034310347465945384

names = []
names.append("H1")
names.append("He4")
names.append("O16")
names.append("O17")
names.append("O18")
names.append("F17")
names.append("F18")
names.append("F19")
names.append("Ne20")
names.append("Ne21")
names.append("Ne22")
names.append("Na21")
names.append("Na22")
names.append("Na23")

def to_composition(Y):
    """Convert an array of molar fractions to a Composition object."""
    from pynucastro import Composition, Nucleus
    nuclei = [Nucleus.from_cache(name) for name in names]
    comp = Composition(nuclei)
    for i, nuc in enumerate(nuclei):
        comp.X[nuc] = Y[i] * A[i]
    return comp


def energy_release(dY):
    """return the energy release in erg/g (/s if dY is actually dY/dt)"""
    enuc = 0.0
    for i, y in enumerate(dY):
        enuc += y * mass[i]
    enuc *= -1*constants.Avogadro
    return enuc

@jitclass([
    ("F17__O17__weak__wc12", numba.float64),
    ("F18__O18__weak__wc12", numba.float64),
    ("Na21__Ne21__weak__wc12", numba.float64),
    ("Na22__Ne22__weak__wc12", numba.float64),
    ("F17__p_O16", numba.float64),
    ("F18__p_O17", numba.float64),
    ("F19__p_O18", numba.float64),
    ("Ne20__p_F19", numba.float64),
    ("Ne20__He4_O16", numba.float64),
    ("Ne21__He4_O17", numba.float64),
    ("Ne22__He4_O18", numba.float64),
    ("Na21__p_Ne20", numba.float64),
    ("Na21__He4_F17", numba.float64),
    ("Na22__p_Ne21", numba.float64),
    ("Na22__He4_F18", numba.float64),
    ("Na23__p_Ne22", numba.float64),
    ("Na23__He4_F19", numba.float64),
    ("p_O16__F17", numba.float64),
    ("He4_O16__Ne20", numba.float64),
    ("p_O17__F18", numba.float64),
    ("He4_O17__Ne21", numba.float64),
    ("p_O18__F19", numba.float64),
    ("He4_O18__Ne22", numba.float64),
    ("He4_F17__Na21", numba.float64),
    ("He4_F18__Na22", numba.float64),
    ("p_F19__Ne20", numba.float64),
    ("He4_F19__Na23", numba.float64),
    ("p_Ne20__Na21", numba.float64),
    ("p_Ne21__Na22", numba.float64),
    ("p_Ne22__Na23", numba.float64),
    ("He4_O16__p_F19", numba.float64),
    ("He4_F17__p_Ne20", numba.float64),
    ("He4_F18__p_Ne21", numba.float64),
    ("p_F19__He4_O16", numba.float64),
    ("He4_F19__p_Ne22", numba.float64),
    ("p_Ne20__He4_F17", numba.float64),
    ("He4_Ne20__p_Na23", numba.float64),
    ("p_Ne21__He4_F18", numba.float64),
    ("p_Ne22__He4_F19", numba.float64),
    ("p_Na23__He4_Ne20", numba.float64),
])
class RateEval:
    def __init__(self):
        self.F17__O17__weak__wc12 = np.nan
        self.F18__O18__weak__wc12 = np.nan
        self.Na21__Ne21__weak__wc12 = np.nan
        self.Na22__Ne22__weak__wc12 = np.nan
        self.F17__p_O16 = np.nan
        self.F18__p_O17 = np.nan
        self.F19__p_O18 = np.nan
        self.Ne20__p_F19 = np.nan
        self.Ne20__He4_O16 = np.nan
        self.Ne21__He4_O17 = np.nan
        self.Ne22__He4_O18 = np.nan
        self.Na21__p_Ne20 = np.nan
        self.Na21__He4_F17 = np.nan
        self.Na22__p_Ne21 = np.nan
        self.Na22__He4_F18 = np.nan
        self.Na23__p_Ne22 = np.nan
        self.Na23__He4_F19 = np.nan
        self.p_O16__F17 = np.nan
        self.He4_O16__Ne20 = np.nan
        self.p_O17__F18 = np.nan
        self.He4_O17__Ne21 = np.nan
        self.p_O18__F19 = np.nan
        self.He4_O18__Ne22 = np.nan
        self.He4_F17__Na21 = np.nan
        self.He4_F18__Na22 = np.nan
        self.p_F19__Ne20 = np.nan
        self.He4_F19__Na23 = np.nan
        self.p_Ne20__Na21 = np.nan
        self.p_Ne21__Na22 = np.nan
        self.p_Ne22__Na23 = np.nan
        self.He4_O16__p_F19 = np.nan
        self.He4_F17__p_Ne20 = np.nan
        self.He4_F18__p_Ne21 = np.nan
        self.p_F19__He4_O16 = np.nan
        self.He4_F19__p_Ne22 = np.nan
        self.p_Ne20__He4_F17 = np.nan
        self.He4_Ne20__p_Na23 = np.nan
        self.p_Ne21__He4_F18 = np.nan
        self.p_Ne22__He4_F19 = np.nan
        self.p_Na23__He4_Ne20 = np.nan

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def F17__O17__weak__wc12(rate_eval, tf):
    # F17 --> O17
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.53318)

    rate_eval.F17__O17__weak__wc12 = rate

@numba.njit()
def F18__O18__weak__wc12(rate_eval, tf):
    # F18 --> O18
    rate = 0.0

    # wc12w
    rate += np.exp(  -9.15982)

    rate_eval.F18__O18__weak__wc12 = rate

@numba.njit()
def Na21__Ne21__weak__wc12(rate_eval, tf):
    # Na21 --> Ne21
    rate = 0.0

    # wc12w
    rate += np.exp(  -3.48003)

    rate_eval.Na21__Ne21__weak__wc12 = rate

@numba.njit()
def Na22__Ne22__weak__wc12(rate_eval, tf):
    # Na22 --> Ne22
    rate = 0.0

    # wc12w
    rate += np.exp(  -18.59)

    rate_eval.Na22__Ne22__weak__wc12 = rate

@numba.njit()
def F17__p_O16(rate_eval, tf):
    # F17 --> p + O16
    rate = 0.0

    # ia08n
    rate += np.exp(  40.9135 + -6.96583*tf.T9i + -16.696*tf.T913i + -1.16252*tf.T913
                  + 0.267703*tf.T9 + -0.0338411*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F17__p_O16 = rate

@numba.njit()
def F18__p_O17(rate_eval, tf):
    # F18 --> p + O17
    rate = 0.0

    # il10r
    rate += np.exp(  33.7037 + -71.2889*tf.T9i + 2.31435*tf.T913
                  + -0.302835*tf.T9 + 0.020133*tf.T953)
    # il10r
    rate += np.exp(  11.2362 + -65.8069*tf.T9i)
    # il10n
    rate += np.exp(  40.2061 + -65.0606*tf.T9i + -16.4035*tf.T913i + 4.31885*tf.T913
                  + -0.709921*tf.T9 + -2.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.F18__p_O17 = rate

@numba.njit()
def F19__p_O18(rate_eval, tf):
    # F19 --> p + O18
    rate = 0.0

    # il10n
    rate += np.exp(  42.8485 + -92.7757*tf.T9i + -16.7246*tf.T913i
                  + -3.0*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  30.2003 + -99.501*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953)
    # il10r
    rate += np.exp(  28.008 + -94.4325*tf.T9i)
    # il10r
    rate += np.exp(  -12.0764 + -93.0204*tf.T9i)

    rate_eval.F19__p_O18 = rate

@numba.njit()
def Ne20__p_F19(rate_eval, tf):
    # Ne20 --> p + F19
    rate = 0.0

    # nacrr
    rate += np.exp(  18.691 + -156.781*tf.T9i + 31.6442*tf.T913i + -58.6563*tf.T913
                  + 67.7365*tf.T9 + -22.9721*tf.T953)
    # nacrr
    rate += np.exp(  36.7036 + -150.75*tf.T9i + -11.3832*tf.T913i + 5.47872*tf.T913
                  + -1.07203*tf.T9 + 0.11196*tf.T953)
    # nacrn
    rate += np.exp(  42.6027 + -149.037*tf.T9i + -18.116*tf.T913i + -1.4622*tf.T913
                  + 6.95113*tf.T9 + -2.90366*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ne20__p_F19 = rate

@numba.njit()
def Ne20__He4_O16(rate_eval, tf):
    # Ne20 --> He4 + O16
    rate = 0.0

    # co10r
    rate += np.exp(  34.2658 + -67.6518*tf.T9i + -3.65925*tf.T913
                  + 0.714224*tf.T9 + -0.00107508*tf.T953)
    # co10r
    rate += np.exp(  28.6431 + -65.246*tf.T9i)
    # co10n
    rate += np.exp(  48.6604 + -54.8875*tf.T9i + -39.7262*tf.T913i + -0.210799*tf.T913
                  + 0.442879*tf.T9 + -0.0797753*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ne20__He4_O16 = rate

@numba.njit()
def Ne21__He4_O17(rate_eval, tf):
    # Ne21 --> He4 + O17
    rate = 0.0

    # be13r
    rate += np.exp(  27.3205 + -91.2722*tf.T9i + 2.87641*tf.T913i + -3.54489*tf.T913
                  + -2.11222e-08*tf.T9 + -3.90649e-09*tf.T953 + 6.25778*tf.lnT9)
    # be13r
    rate += np.exp(  0.0906657 + -90.782*tf.T9i + 123.363*tf.T913i + -87.4351*tf.T913
                  + -3.40974e-06*tf.T9 + -57.0469*tf.T953 + 83.7218*tf.lnT9)
    # be13r
    rate += np.exp(  -91.954 + -98.9487*tf.T9i + 3.31162e-08*tf.T913i + 130.258*tf.T913
                  + -7.92551e-05*tf.T9 + -4.13772*tf.T953 + -41.2753*tf.lnT9)

    rate_eval.Ne21__He4_O17 = rate

@numba.njit()
def Ne22__He4_O18(rate_eval, tf):
    # Ne22 --> He4 + O18
    rate = 0.0

    # il10r
    rate += np.exp(  39.7659 + -143.24*tf.T9i)
    # il10r
    rate += np.exp(  106.996 + -113.779*tf.T9i + -44.3823*tf.T913i + -46.6617*tf.T913
                  + 7.88059*tf.T9 + -0.590829*tf.T953)
    # il10r
    rate += np.exp(  -7.12154 + -114.197*tf.T9i)
    # il10r
    rate += np.exp(  -56.5125 + -112.87*tf.T9i)

    rate_eval.Ne22__He4_O18 = rate

@numba.njit()
def Na21__p_Ne20(rate_eval, tf):
    # Na21 --> p + Ne20
    rate = 0.0

    # ly18 
    rate += np.exp(  195320.0 + -89.3596*tf.T9i + 21894.7*tf.T913i + -319153.0*tf.T913
                  + 224369.0*tf.T9 + -188049.0*tf.T953 + 48704.9*tf.lnT9)
    # ly18 
    rate += np.exp(  230.123 + -28.3722*tf.T9i + 15.325*tf.T913i + -294.859*tf.T913
                  + 107.692*tf.T9 + -46.2072*tf.T953 + 59.3398*tf.lnT9)
    # ly18 
    rate += np.exp(  28.0772 + -37.0575*tf.T9i + 20.5893*tf.T913i + -17.5841*tf.T913
                  + 0.243226*tf.T9 + -0.000231418*tf.T953 + 14.3398*tf.lnT9)
    # ly18 
    rate += np.exp(  252.265 + -32.6731*tf.T9i + 258.57*tf.T913i + -506.387*tf.T913
                  + 22.1576*tf.T9 + -0.721182*tf.T953 + 231.788*tf.lnT9)

    rate_eval.Na21__p_Ne20 = rate

@numba.njit()
def Na21__He4_F17(rate_eval, tf):
    # Na21 --> He4 + F17
    rate = 0.0

    # rpsmr
    rate += np.exp(  66.3334 + -77.8653*tf.T9i + 15.559*tf.T913i + -68.3231*tf.T913
                  + 2.54275*tf.T9 + -0.0989207*tf.T953 + 38.3877*tf.lnT9)

    rate_eval.Na21__He4_F17 = rate

@numba.njit()
def Na22__p_Ne21(rate_eval, tf):
    # Na22 --> p + Ne21
    rate = 0.0

    # il10r
    rate += np.exp(  -16.4098 + -82.4235*tf.T9i + 21.1176*tf.T913i + 34.0411*tf.T913
                  + -4.45593*tf.T9 + 0.328613*tf.T953)
    # il10r
    rate += np.exp(  24.8334 + -79.6093*tf.T9i)
    # il10r
    rate += np.exp(  -24.579 + -78.4059*tf.T9i)
    # il10n
    rate += np.exp(  42.146 + -78.2097*tf.T9i + -19.2096*tf.T913i
                  + -1.0*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Na22__p_Ne21 = rate

@numba.njit()
def Na22__He4_F18(rate_eval, tf):
    # Na22 --> He4 + F18
    rate = 0.0

    # rpsmr
    rate += np.exp(  59.3224 + -100.236*tf.T9i + 18.8956*tf.T913i + -65.6134*tf.T913
                  + 1.71114*tf.T9 + -0.0260999*tf.T953 + 39.3396*tf.lnT9)

    rate_eval.Na22__He4_F18 = rate

@numba.njit()
def Na23__p_Ne22(rate_eval, tf):
    # Na23 --> p + Ne22
    rate = 0.0

    # ke17r
    rate += np.exp(  18.2467 + -104.673*tf.T9i
                  + -2.79964*tf.lnT9)
    # ke17r
    rate += np.exp(  21.6534 + -103.776*tf.T9i
                  + 1.18923*tf.lnT9)
    # ke17r
    rate += np.exp(  0.818178 + -102.466*tf.T9i
                  + 0.009812*tf.lnT9)
    # ke17r
    rate += np.exp(  18.1624 + -102.855*tf.T9i
                  + 4.73558*tf.lnT9)
    # ke17r
    rate += np.exp(  36.29 + -110.779*tf.T9i
                  + 0.732533*tf.lnT9)
    # ke17r
    rate += np.exp(  33.8935 + -106.655*tf.T9i
                  + 1.65623*tf.lnT9)

    rate_eval.Na23__p_Ne22 = rate

@numba.njit()
def Na23__He4_F19(rate_eval, tf):
    # Na23 --> He4 + F19
    rate = 0.0

    # rpsmr
    rate += np.exp(  76.8979 + -123.578*tf.T9i + 39.7219*tf.T913i + -100.401*tf.T913
                  + 3.15808*tf.T9 + -0.0629822*tf.T953 + 55.9823*tf.lnT9)

    rate_eval.Na23__He4_F19 = rate

@numba.njit()
def p_O16__F17(rate_eval, tf):
    # O16 + p --> F17
    rate = 0.0

    # ia08n
    rate += np.exp(  19.0904 + -16.696*tf.T913i + -1.16252*tf.T913
                  + 0.267703*tf.T9 + -0.0338411*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_O16__F17 = rate

@numba.njit()
def He4_O16__Ne20(rate_eval, tf):
    # O16 + He4 --> Ne20
    rate = 0.0

    # co10r
    rate += np.exp(  9.50848 + -12.7643*tf.T9i + -3.65925*tf.T913
                  + 0.714224*tf.T9 + -0.00107508*tf.T953 + -1.5*tf.lnT9)
    # co10r
    rate += np.exp(  3.88571 + -10.3585*tf.T9i
                  + -1.5*tf.lnT9)
    # co10n
    rate += np.exp(  23.903 + -39.7262*tf.T913i + -0.210799*tf.T913
                  + 0.442879*tf.T9 + -0.0797753*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O16__Ne20 = rate

@numba.njit()
def p_O17__F18(rate_eval, tf):
    # O17 + p --> F18
    rate = 0.0

    # il10n
    rate += np.exp(  15.8929 + -16.4035*tf.T913i + 4.31885*tf.T913
                  + -0.709921*tf.T9 + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  9.39048 + -6.22828*tf.T9i + 2.31435*tf.T913
                  + -0.302835*tf.T9 + 0.020133*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -13.077 + -0.746296*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O17__F18 = rate

@numba.njit()
def He4_O17__Ne21(rate_eval, tf):
    # O17 + He4 --> Ne21
    rate = 0.0

    # be13r
    rate += np.exp(  -25.0898 + -5.50926*tf.T9i + 123.363*tf.T913i + -87.4351*tf.T913
                  + -3.40974e-06*tf.T9 + -57.0469*tf.T953 + 82.2218*tf.lnT9)
    # be13r
    rate += np.exp(  -117.134 + -13.6759*tf.T9i + 3.31162e-08*tf.T913i + 130.258*tf.T913
                  + -7.92551e-05*tf.T9 + -4.13772*tf.T953 + -42.7753*tf.lnT9)
    # be13r
    rate += np.exp(  2.14 + -5.99952*tf.T9i + 2.87641*tf.T913i + -3.54489*tf.T913
                  + -2.11222e-08*tf.T9 + -3.90649e-09*tf.T953 + 4.75778*tf.lnT9)

    rate_eval.He4_O17__Ne21 = rate

@numba.njit()
def p_O18__F19(rate_eval, tf):
    # O18 + p --> F19
    rate = 0.0

    # il10r
    rate += np.exp(  -35.0079 + -0.244743*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.917 + -16.7246*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  7.26876 + -6.7253*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.07648 + -1.65681*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O18__F19 = rate

@numba.njit()
def He4_O18__Ne22(rate_eval, tf):
    # O18 + He4 --> Ne22
    rate = 0.0

    # il10r
    rate += np.exp(  -81.3036 + -0.676112*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  14.9748 + -31.0468*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  82.2053 + -1.58534*tf.T9i + -44.3823*tf.T913i + -46.6617*tf.T913
                  + 7.88059*tf.T9 + -0.590829*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -31.9126 + -2.00306*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O18__Ne22 = rate

@numba.njit()
def He4_F17__Na21(rate_eval, tf):
    # F17 + He4 --> Na21
    rate = 0.0

    # rpsmr
    rate += np.exp(  41.1529 + -1.72817*tf.T9i + 15.559*tf.T913i + -68.3231*tf.T913
                  + 2.54275*tf.T9 + -0.0989207*tf.T953 + 36.8877*tf.lnT9)

    rate_eval.He4_F17__Na21 = rate

@numba.njit()
def He4_F18__Na22(rate_eval, tf):
    # F18 + He4 --> Na22
    rate = 0.0

    # rpsmr
    rate += np.exp(  35.3786 + -1.82957*tf.T9i + 18.8956*tf.T913i + -65.6134*tf.T913
                  + 1.71114*tf.T9 + -0.0260999*tf.T953 + 37.8396*tf.lnT9)

    rate_eval.He4_F18__Na22 = rate

@numba.njit()
def p_F19__Ne20(rate_eval, tf):
    # F19 + p --> Ne20
    rate = 0.0

    # nacrr
    rate += np.exp(  -5.63093 + -7.74414*tf.T9i + 31.6442*tf.T913i + -58.6563*tf.T913
                  + 67.7365*tf.T9 + -22.9721*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  12.3816 + -1.71383*tf.T9i + -11.3832*tf.T913i + 5.47872*tf.T913
                  + -1.07203*tf.T9 + 0.11196*tf.T953 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.2807 + -18.116*tf.T913i + -1.4622*tf.T913
                  + 6.95113*tf.T9 + -2.90366*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_F19__Ne20 = rate

@numba.njit()
def He4_F19__Na23(rate_eval, tf):
    # F19 + He4 --> Na23
    rate = 0.0

    # rpsmr
    rate += np.exp(  52.7856 + -2.11408*tf.T9i + 39.7219*tf.T913i + -100.401*tf.T913
                  + 3.15808*tf.T9 + -0.0629822*tf.T953 + 54.4823*tf.lnT9)

    rate_eval.He4_F19__Na23 = rate

@numba.njit()
def p_Ne20__Na21(rate_eval, tf):
    # Ne20 + p --> Na21
    rate = 0.0

    # ly18 
    rate += np.exp(  230.019 + -4.45358*tf.T9i + 258.57*tf.T913i + -506.387*tf.T913
                  + 22.1576*tf.T9 + -0.721182*tf.T953 + 230.288*tf.lnT9)
    # ly18 
    rate += np.exp(  195297.0 + -61.14*tf.T9i + 21894.7*tf.T913i + -319153.0*tf.T913
                  + 224369.0*tf.T9 + -188049.0*tf.T953 + 48703.4*tf.lnT9)
    # ly18 
    rate += np.exp(  207.877 + -0.152711*tf.T9i + 15.325*tf.T913i + -294.859*tf.T913
                  + 107.692*tf.T9 + -46.2072*tf.T953 + 57.8398*tf.lnT9)
    # ly18 
    rate += np.exp(  5.83103 + -8.838*tf.T9i + 20.5893*tf.T913i + -17.5841*tf.T913
                  + 0.243226*tf.T9 + -0.000231418*tf.T953 + 12.8398*tf.lnT9)

    rate_eval.p_Ne20__Na21 = rate

@numba.njit()
def p_Ne21__Na22(rate_eval, tf):
    # Ne21 + p --> Na22
    rate = 0.0

    # il10r
    rate += np.exp(  -47.6554 + -0.19618*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.0696 + -19.2096*tf.T913i
                  + -1.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  -39.4862 + -4.21385*tf.T9i + 21.1176*tf.T913i + 34.0411*tf.T913
                  + -4.45593*tf.T9 + 0.328613*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  1.75704 + -1.39957*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Ne21__Na22 = rate

@numba.njit()
def p_Ne22__Na23(rate_eval, tf):
    # Ne22 + p --> Na23
    rate = 0.0

    # ke17r
    rate += np.exp(  -4.00597 + -2.6179*tf.T9i
                  + -4.29964*tf.lnT9)
    # ke17r
    rate += np.exp(  -0.599331 + -1.72007*tf.T9i
                  + -0.310765*tf.lnT9)
    # ke17r
    rate += np.exp(  -21.4345 + -0.410962*tf.T9i
                  + -1.49019*tf.lnT9)
    # ke17r
    rate += np.exp(  -4.09035 + -0.799756*tf.T9i
                  + 3.23558*tf.lnT9)
    # ke17r
    rate += np.exp(  14.0373 + -8.72377*tf.T9i
                  + -0.767467*tf.lnT9)
    # ke17r
    rate += np.exp(  11.6408 + -4.59936*tf.T9i
                  + 0.156226*tf.lnT9)

    rate_eval.p_Ne22__Na23 = rate

@numba.njit()
def He4_O16__p_F19(rate_eval, tf):
    # O16 + He4 --> p + F19
    rate = 0.0

    # nacr 
    rate += np.exp(  -53.1397 + -94.2866*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  25.8562 + -94.1589*tf.T9i + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9232 + -97.4449*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  14.7601 + -97.9108*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  7.80363 + -96.6272*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O16__p_F19 = rate

@numba.njit()
def He4_F17__p_Ne20(rate_eval, tf):
    # F17 + He4 --> p + Ne20
    rate = 0.0

    # nacr 
    rate += np.exp(  38.6287 + -43.18*tf.T913i + 4.46827*tf.T913
                  + -1.63915*tf.T9 + 0.123483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_F17__p_Ne20 = rate

@numba.njit()
def He4_F18__p_Ne21(rate_eval, tf):
    # F18 + He4 --> p + Ne21
    rate = 0.0

    # rpsmr
    rate += np.exp(  49.7863 + -1.84559*tf.T9i + 21.4461*tf.T913i + -73.252*tf.T913
                  + 2.42329*tf.T9 + -0.077278*tf.T953 + 40.7604*tf.lnT9)

    rate_eval.He4_F18__p_Ne21 = rate

@numba.njit()
def p_F19__He4_O16(rate_eval, tf):
    # F19 + p --> He4 + O16
    rate = 0.0

    # nacr 
    rate += np.exp(  8.239 + -2.46828*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  -52.7043 + -0.12765*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  26.2916 + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  14.3586 + -3.286*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  15.1955 + -3.75185*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_F19__He4_O16 = rate

@numba.njit()
def He4_F19__p_Ne22(rate_eval, tf):
    # F19 + He4 --> p + Ne22
    rate = 0.0

    # da18r
    rate += np.exp(  29430.6 + -133.026*tf.T9i + 12625.1*tf.T913i + -49107.1*tf.T913
                  + 9227.53*tf.T9 + -2086.65*tf.T953 + 14520.2*tf.lnT9)
    # da18r
    rate += np.exp(  52.9317 + -2.8444*tf.T9i + -38.7722*tf.T913i + -13.3654*tf.T913
                  + 0.863648*tf.T9 + -0.0451491*tf.T953 + 1.33333*tf.lnT9)
    # da18r
    rate += np.exp(  51.6709 + -45.7808*tf.T9i + -34.5008*tf.T913i + 56.9316*tf.T913
                  + 2.09613*tf.T9 + -32.496*tf.T953 + 0.333333*tf.lnT9)

    rate_eval.He4_F19__p_Ne22 = rate

@numba.njit()
def p_Ne20__He4_F17(rate_eval, tf):
    # Ne20 + p --> He4 + F17
    rate = 0.0

    # nacr 
    rate += np.exp(  41.563 + -47.9266*tf.T9i + -43.18*tf.T913i + 4.46827*tf.T913
                  + -1.63915*tf.T9 + 0.123483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ne20__He4_F17 = rate

@numba.njit()
def He4_Ne20__p_Na23(rate_eval, tf):
    # Ne20 + He4 --> p + Na23
    rate = 0.0

    # il10r
    rate += np.exp(  0.227472 + -29.4348*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.1852 + -27.5738*tf.T9i + -20.0024*tf.T913i + 11.5988*tf.T913
                  + -1.37398*tf.T9 + -1.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  -6.37772 + -29.8896*tf.T9i + 19.7297*tf.T913
                  + -2.20987*tf.T9 + 0.153374*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_Ne20__p_Na23 = rate

@numba.njit()
def p_Ne21__He4_F18(rate_eval, tf):
    # Ne21 + p --> He4 + F18
    rate = 0.0

    # rpsmr
    rate += np.exp(  50.6536 + -22.049*tf.T9i + 21.4461*tf.T913i + -73.252*tf.T913
                  + 2.42329*tf.T9 + -0.077278*tf.T953 + 40.7604*tf.lnT9)

    rate_eval.p_Ne21__He4_F18 = rate

@numba.njit()
def p_Ne22__He4_F19(rate_eval, tf):
    # Ne22 + p --> He4 + F19
    rate = 0.0

    # da18r
    rate += np.exp(  53.5304 + -65.1991*tf.T9i + -34.5008*tf.T913i + 56.9316*tf.T913
                  + 2.09613*tf.T9 + -32.496*tf.T953 + 0.333333*tf.lnT9)
    # da18r
    rate += np.exp(  29432.5 + -152.444*tf.T9i + 12625.1*tf.T913i + -49107.1*tf.T913
                  + 9227.53*tf.T9 + -2086.65*tf.T953 + 14520.2*tf.lnT9)
    # da18r
    rate += np.exp(  54.7912 + -22.2627*tf.T9i + -38.7722*tf.T913i + -13.3654*tf.T913
                  + 0.863648*tf.T9 + -0.0451491*tf.T953 + 1.33333*tf.lnT9)

    rate_eval.p_Ne22__He4_F19 = rate

@numba.njit()
def p_Na23__He4_Ne20(rate_eval, tf):
    # Na23 + p --> He4 + Ne20
    rate = 0.0

    # il10r
    rate += np.exp(  -6.58736 + -2.31577*tf.T9i + 19.7297*tf.T913
                  + -2.20987*tf.T9 + 0.153374*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  0.0178295 + -1.86103*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  18.9756 + -20.0024*tf.T913i + 11.5988*tf.T913
                  + -1.37398*tf.T9 + -1.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Na23__He4_Ne20 = rate

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    F17__O17__weak__wc12(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    Na21__Ne21__weak__wc12(rate_eval, tf)
    Na22__Ne22__weak__wc12(rate_eval, tf)
    F17__p_O16(rate_eval, tf)
    F18__p_O17(rate_eval, tf)
    F19__p_O18(rate_eval, tf)
    Ne20__p_F19(rate_eval, tf)
    Ne20__He4_O16(rate_eval, tf)
    Ne21__He4_O17(rate_eval, tf)
    Ne22__He4_O18(rate_eval, tf)
    Na21__p_Ne20(rate_eval, tf)
    Na21__He4_F17(rate_eval, tf)
    Na22__p_Ne21(rate_eval, tf)
    Na22__He4_F18(rate_eval, tf)
    Na23__p_Ne22(rate_eval, tf)
    Na23__He4_F19(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    He4_O16__Ne20(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    He4_O17__Ne21(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    He4_O18__Ne22(rate_eval, tf)
    He4_F17__Na21(rate_eval, tf)
    He4_F18__Na22(rate_eval, tf)
    p_F19__Ne20(rate_eval, tf)
    He4_F19__Na23(rate_eval, tf)
    p_Ne20__Na21(rate_eval, tf)
    p_Ne21__Na22(rate_eval, tf)
    p_Ne22__Na23(rate_eval, tf)
    He4_O16__p_F19(rate_eval, tf)
    He4_F17__p_Ne20(rate_eval, tf)
    He4_F18__p_Ne21(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    He4_F19__p_Ne22(rate_eval, tf)
    p_Ne20__He4_F17(rate_eval, tf)
    He4_Ne20__p_Na23(rate_eval, tf)
    p_Ne21__He4_F18(rate_eval, tf)
    p_Ne22__He4_F19(rate_eval, tf)
    p_Na23__He4_Ne20(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__Ne20 *= scor
        rate_eval.He4_O16__p_F19 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__F18 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O17__Ne21 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O18__Ne22 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F17__Na21 *= scor
        rate_eval.He4_F17__p_Ne20 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F18__Na22 *= scor
        rate_eval.He4_F18__p_Ne21 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__Ne20 *= scor
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F19__Na23 *= scor
        rate_eval.He4_F19__p_Ne22 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne20__Na21 *= scor
        rate_eval.p_Ne20__He4_F17 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne21__Na22 *= scor
        rate_eval.p_Ne21__He4_F18 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne22__Na23 *= scor
        rate_eval.p_Ne22__He4_F19 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne20__p_Na23 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na23__He4_Ne20 *= scor

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jp] = (
       -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__Na21
       -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__Na22
       -rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__Na23
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18
       -rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__He4_F19
       -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       +Y[jf17]*rate_eval.F17__p_O16
       +Y[jf18]*rate_eval.F18__p_O17
       +Y[jf19]*rate_eval.F19__p_O18
       +Y[jne20]*rate_eval.Ne20__p_F19
       +Y[jna21]*rate_eval.Na21__p_Ne20
       +Y[jna22]*rate_eval.Na22__p_Ne21
       +Y[jna23]*rate_eval.Na23__p_Ne22
       +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21
       +rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__p_Ne22
       +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       )

    dYdt[jhe4] = (
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__Ne21
       -rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__Ne22
       -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__Na21
       -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__Na22
       -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__Na23
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20
       -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21
       -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__p_Ne22
       -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       +Y[jne20]*rate_eval.Ne20__He4_O16
       +Y[jne21]*rate_eval.Ne21__He4_O17
       +Y[jne22]*rate_eval.Ne22__He4_O18
       +Y[jna21]*rate_eval.Na21__He4_F17
       +Y[jna22]*rate_eval.Na22__He4_F18
       +Y[jna23]*rate_eval.Na23__He4_F19
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18
       +rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__He4_F19
       +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       )

    dYdt[jo16] = (
       -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       +Y[jf17]*rate_eval.F17__p_O16
       +Y[jne20]*rate_eval.Ne20__He4_O16
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    dYdt[jo17] = (
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__Ne21
       +Y[jf17]*rate_eval.F17__O17__weak__wc12
       +Y[jf18]*rate_eval.F18__p_O17
       +Y[jne21]*rate_eval.Ne21__He4_O17
       )

    dYdt[jo18] = (
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__Ne22
       +Y[jf18]*rate_eval.F18__O18__weak__wc12
       +Y[jf19]*rate_eval.F19__p_O18
       +Y[jne22]*rate_eval.Ne22__He4_O18
       )

    dYdt[jf17] = (
       -Y[jf17]*rate_eval.F17__O17__weak__wc12
       -Y[jf17]*rate_eval.F17__p_O16
       -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__Na21
       -rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +Y[jna21]*rate_eval.Na21__He4_F17
       +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    dYdt[jf18] = (
       -Y[jf18]*rate_eval.F18__O18__weak__wc12
       -Y[jf18]*rate_eval.F18__p_O17
       -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__Na22
       -rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21
       +Y[jna22]*rate_eval.Na22__He4_F18
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18
       )

    dYdt[jf19] = (
       -Y[jf19]*rate_eval.F19__p_O18
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__Na23
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__p_Ne22
       +Y[jne20]*rate_eval.Ne20__p_F19
       +Y[jna23]*rate_eval.Na23__He4_F19
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__He4_F19
       )

    dYdt[jne20] = (
       -Y[jne20]*rate_eval.Ne20__p_F19
       -Y[jne20]*rate_eval.Ne20__He4_O16
       -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__Na21
       -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       +Y[jna21]*rate_eval.Na21__p_Ne20
       +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20
       +rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       )

    dYdt[jne21] = (
       -Y[jne21]*rate_eval.Ne21__He4_O17
       -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__Na22
       -rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__He4_F18
       +Y[jna21]*rate_eval.Na21__Ne21__weak__wc12
       +Y[jna22]*rate_eval.Na22__p_Ne21
       +rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__Ne21
       +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__p_Ne21
       )

    dYdt[jne22] = (
       -Y[jne22]*rate_eval.Ne22__He4_O18
       -rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__Na23
       -rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__He4_F19
       +Y[jna22]*rate_eval.Na22__Ne22__weak__wc12
       +Y[jna23]*rate_eval.Na23__p_Ne22
       +rho*Y[jhe4]*Y[jo18]*rate_eval.He4_O18__Ne22
       +rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__p_Ne22
       )

    dYdt[jna21] = (
       -Y[jna21]*rate_eval.Na21__Ne21__weak__wc12
       -Y[jna21]*rate_eval.Na21__p_Ne20
       -Y[jna21]*rate_eval.Na21__He4_F17
       +rho*Y[jhe4]*Y[jf17]*rate_eval.He4_F17__Na21
       +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__Na21
       )

    dYdt[jna22] = (
       -Y[jna22]*rate_eval.Na22__Ne22__weak__wc12
       -Y[jna22]*rate_eval.Na22__p_Ne21
       -Y[jna22]*rate_eval.Na22__He4_F18
       +rho*Y[jhe4]*Y[jf18]*rate_eval.He4_F18__Na22
       +rho*Y[jp]*Y[jne21]*rate_eval.p_Ne21__Na22
       )

    dYdt[jna23] = (
       -Y[jna23]*rate_eval.Na23__p_Ne22
       -Y[jna23]*rate_eval.Na23__He4_F19
       -rho*Y[jp]*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       +rho*Y[jhe4]*Y[jf19]*rate_eval.He4_F19__Na23
       +rho*Y[jp]*Y[jne22]*rate_eval.p_Ne22__Na23
       +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    F17__O17__weak__wc12(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    Na21__Ne21__weak__wc12(rate_eval, tf)
    Na22__Ne22__weak__wc12(rate_eval, tf)
    F17__p_O16(rate_eval, tf)
    F18__p_O17(rate_eval, tf)
    F19__p_O18(rate_eval, tf)
    Ne20__p_F19(rate_eval, tf)
    Ne20__He4_O16(rate_eval, tf)
    Ne21__He4_O17(rate_eval, tf)
    Ne22__He4_O18(rate_eval, tf)
    Na21__p_Ne20(rate_eval, tf)
    Na21__He4_F17(rate_eval, tf)
    Na22__p_Ne21(rate_eval, tf)
    Na22__He4_F18(rate_eval, tf)
    Na23__p_Ne22(rate_eval, tf)
    Na23__He4_F19(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    He4_O16__Ne20(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    He4_O17__Ne21(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    He4_O18__Ne22(rate_eval, tf)
    He4_F17__Na21(rate_eval, tf)
    He4_F18__Na22(rate_eval, tf)
    p_F19__Ne20(rate_eval, tf)
    He4_F19__Na23(rate_eval, tf)
    p_Ne20__Na21(rate_eval, tf)
    p_Ne21__Na22(rate_eval, tf)
    p_Ne22__Na23(rate_eval, tf)
    He4_O16__p_F19(rate_eval, tf)
    He4_F17__p_Ne20(rate_eval, tf)
    He4_F18__p_Ne21(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    He4_F19__p_Ne22(rate_eval, tf)
    p_Ne20__He4_F17(rate_eval, tf)
    He4_Ne20__p_Na23(rate_eval, tf)
    p_Ne21__He4_F18(rate_eval, tf)
    p_Ne22__He4_F19(rate_eval, tf)
    p_Na23__He4_Ne20(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__Ne20 *= scor
        rate_eval.He4_O16__p_F19 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__F18 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O17__Ne21 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O18__Ne22 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F17__Na21 *= scor
        rate_eval.He4_F17__p_Ne20 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F18__Na22 *= scor
        rate_eval.He4_F18__p_Ne21 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__Ne20 *= scor
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_F19__Na23 *= scor
        rate_eval.He4_F19__p_Ne22 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne20__Na21 *= scor
        rate_eval.p_Ne20__He4_F17 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 21)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne21__Na22 *= scor
        rate_eval.p_Ne21__He4_F18 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 22)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne22__Na23 *= scor
        rate_eval.p_Ne22__He4_F19 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne20__p_Na23 *= scor

        scn_fac = ScreenFactors(1, 1, 11, 23)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Na23__He4_Ne20 *= scor

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jp, jp] = (
       -rho*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jne20]*rate_eval.p_Ne20__Na21
       -rho*Y[jne21]*rate_eval.p_Ne21__Na22
       -rho*Y[jne22]*rate_eval.p_Ne22__Na23
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       -rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       -rho*Y[jne22]*rate_eval.p_Ne22__He4_F19
       -rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jp, jhe4] = (
       +rho*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       +rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       +rho*Y[jf19]*rate_eval.He4_F19__p_Ne22
       +rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jp, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       +rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jp, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jp, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jp, jf17] = (
       +rate_eval.F17__p_O16
       +rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jp, jf18] = (
       +rate_eval.F18__p_O17
       +rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jp, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       +rate_eval.F19__p_O18
       +rho*Y[jhe4]*rate_eval.He4_F19__p_Ne22
       )

    jac[jp, jne20] = (
       -rho*Y[jp]*rate_eval.p_Ne20__Na21
       -rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       +rate_eval.Ne20__p_F19
       +rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jp, jne21] = (
       -rho*Y[jp]*rate_eval.p_Ne21__Na22
       -rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       )

    jac[jp, jne22] = (
       -rho*Y[jp]*rate_eval.p_Ne22__Na23
       -rho*Y[jp]*rate_eval.p_Ne22__He4_F19
       )

    jac[jp, jna21] = (
       +rate_eval.Na21__p_Ne20
       )

    jac[jp, jna22] = (
       +rate_eval.Na22__p_Ne21
       )

    jac[jp, jna23] = (
       -rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       +rate_eval.Na23__p_Ne22
       )

    jac[jhe4, jp] = (
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       +rho*Y[jne22]*rate_eval.p_Ne22__He4_F19
       +rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jhe4, jhe4] = (
       -rho*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jo17]*rate_eval.He4_O17__Ne21
       -rho*Y[jo18]*rate_eval.He4_O18__Ne22
       -rho*Y[jf17]*rate_eval.He4_F17__Na21
       -rho*Y[jf18]*rate_eval.He4_F18__Na22
       -rho*Y[jf19]*rate_eval.He4_F19__Na23
       -rho*Y[jo16]*rate_eval.He4_O16__p_F19
       -rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       -rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       -rho*Y[jf19]*rate_eval.He4_F19__p_Ne22
       -rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jhe4, jo16] = (
       -rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       -rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jhe4, jo17] = (
       -rho*Y[jhe4]*rate_eval.He4_O17__Ne21
       )

    jac[jhe4, jo18] = (
       -rho*Y[jhe4]*rate_eval.He4_O18__Ne22
       )

    jac[jhe4, jf17] = (
       -rho*Y[jhe4]*rate_eval.He4_F17__Na21
       -rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jhe4, jf18] = (
       -rho*Y[jhe4]*rate_eval.He4_F18__Na22
       -rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jhe4, jf19] = (
       -rho*Y[jhe4]*rate_eval.He4_F19__Na23
       -rho*Y[jhe4]*rate_eval.He4_F19__p_Ne22
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jhe4, jne20] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       +rate_eval.Ne20__He4_O16
       +rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jhe4, jne21] = (
       +rate_eval.Ne21__He4_O17
       +rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       )

    jac[jhe4, jne22] = (
       +rate_eval.Ne22__He4_O18
       +rho*Y[jp]*rate_eval.p_Ne22__He4_F19
       )

    jac[jhe4, jna21] = (
       +rate_eval.Na21__He4_F17
       )

    jac[jhe4, jna22] = (
       +rate_eval.Na22__He4_F18
       )

    jac[jhe4, jna23] = (
       +rate_eval.Na23__He4_F19
       +rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jo16, jp] = (
       -rho*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    jac[jo16, jhe4] = (
       -rho*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jo16]*rate_eval.He4_O16__p_F19
       )

    jac[jo16, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       -rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       -rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jo16, jf17] = (
       +rate_eval.F17__p_O16
       )

    jac[jo16, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jo16, jne20] = (
       +rate_eval.Ne20__He4_O16
       )

    jac[jo17, jp] = (
       -rho*Y[jo17]*rate_eval.p_O17__F18
       )

    jac[jo17, jhe4] = (
       -rho*Y[jo17]*rate_eval.He4_O17__Ne21
       )

    jac[jo17, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__F18
       -rho*Y[jhe4]*rate_eval.He4_O17__Ne21
       )

    jac[jo17, jf17] = (
       +rate_eval.F17__O17__weak__wc12
       )

    jac[jo17, jf18] = (
       +rate_eval.F18__p_O17
       )

    jac[jo17, jne21] = (
       +rate_eval.Ne21__He4_O17
       )

    jac[jo18, jp] = (
       -rho*Y[jo18]*rate_eval.p_O18__F19
       )

    jac[jo18, jhe4] = (
       -rho*Y[jo18]*rate_eval.He4_O18__Ne22
       )

    jac[jo18, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jhe4]*rate_eval.He4_O18__Ne22
       )

    jac[jo18, jf18] = (
       +rate_eval.F18__O18__weak__wc12
       )

    jac[jo18, jf19] = (
       +rate_eval.F19__p_O18
       )

    jac[jo18, jne22] = (
       +rate_eval.Ne22__He4_O18
       )

    jac[jf17, jp] = (
       +rho*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    jac[jf17, jhe4] = (
       -rho*Y[jf17]*rate_eval.He4_F17__Na21
       -rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       )

    jac[jf17, jo16] = (
       +rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jf17, jf17] = (
       -rate_eval.F17__O17__weak__wc12
       -rate_eval.F17__p_O16
       -rho*Y[jhe4]*rate_eval.He4_F17__Na21
       -rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jf17, jne20] = (
       +rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jf17, jna21] = (
       +rate_eval.Na21__He4_F17
       )

    jac[jf18, jp] = (
       +rho*Y[jo17]*rate_eval.p_O17__F18
       +rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       )

    jac[jf18, jhe4] = (
       -rho*Y[jf18]*rate_eval.He4_F18__Na22
       -rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       )

    jac[jf18, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jf18, jf18] = (
       -rate_eval.F18__O18__weak__wc12
       -rate_eval.F18__p_O17
       -rho*Y[jhe4]*rate_eval.He4_F18__Na22
       -rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jf18, jne21] = (
       +rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       )

    jac[jf18, jna22] = (
       +rate_eval.Na22__He4_F18
       )

    jac[jf19, jp] = (
       -rho*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jo18]*rate_eval.p_O18__F19
       +rho*Y[jne22]*rate_eval.p_Ne22__He4_F19
       )

    jac[jf19, jhe4] = (
       -rho*Y[jf19]*rate_eval.He4_F19__Na23
       -rho*Y[jf19]*rate_eval.He4_F19__p_Ne22
       +rho*Y[jo16]*rate_eval.He4_O16__p_F19
       )

    jac[jf19, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jf19, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jf19, jf19] = (
       -rate_eval.F19__p_O18
       -rho*Y[jp]*rate_eval.p_F19__Ne20
       -rho*Y[jhe4]*rate_eval.He4_F19__Na23
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       -rho*Y[jhe4]*rate_eval.He4_F19__p_Ne22
       )

    jac[jf19, jne20] = (
       +rate_eval.Ne20__p_F19
       )

    jac[jf19, jne22] = (
       +rho*Y[jp]*rate_eval.p_Ne22__He4_F19
       )

    jac[jf19, jna23] = (
       +rate_eval.Na23__He4_F19
       )

    jac[jne20, jp] = (
       -rho*Y[jne20]*rate_eval.p_Ne20__Na21
       -rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jf19]*rate_eval.p_F19__Ne20
       +rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jne20, jhe4] = (
       -rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       +rho*Y[jo16]*rate_eval.He4_O16__Ne20
       +rho*Y[jf17]*rate_eval.He4_F17__p_Ne20
       )

    jac[jne20, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       )

    jac[jne20, jf17] = (
       +rho*Y[jhe4]*rate_eval.He4_F17__p_Ne20
       )

    jac[jne20, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__Ne20
       )

    jac[jne20, jne20] = (
       -rate_eval.Ne20__p_F19
       -rate_eval.Ne20__He4_O16
       -rho*Y[jp]*rate_eval.p_Ne20__Na21
       -rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       -rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jne20, jna21] = (
       +rate_eval.Na21__p_Ne20
       )

    jac[jne20, jna23] = (
       +rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       )

    jac[jne21, jp] = (
       -rho*Y[jne21]*rate_eval.p_Ne21__Na22
       -rho*Y[jne21]*rate_eval.p_Ne21__He4_F18
       )

    jac[jne21, jhe4] = (
       +rho*Y[jo17]*rate_eval.He4_O17__Ne21
       +rho*Y[jf18]*rate_eval.He4_F18__p_Ne21
       )

    jac[jne21, jo17] = (
       +rho*Y[jhe4]*rate_eval.He4_O17__Ne21
       )

    jac[jne21, jf18] = (
       +rho*Y[jhe4]*rate_eval.He4_F18__p_Ne21
       )

    jac[jne21, jne21] = (
       -rate_eval.Ne21__He4_O17
       -rho*Y[jp]*rate_eval.p_Ne21__Na22
       -rho*Y[jp]*rate_eval.p_Ne21__He4_F18
       )

    jac[jne21, jna21] = (
       +rate_eval.Na21__Ne21__weak__wc12
       )

    jac[jne21, jna22] = (
       +rate_eval.Na22__p_Ne21
       )

    jac[jne22, jp] = (
       -rho*Y[jne22]*rate_eval.p_Ne22__Na23
       -rho*Y[jne22]*rate_eval.p_Ne22__He4_F19
       )

    jac[jne22, jhe4] = (
       +rho*Y[jo18]*rate_eval.He4_O18__Ne22
       +rho*Y[jf19]*rate_eval.He4_F19__p_Ne22
       )

    jac[jne22, jo18] = (
       +rho*Y[jhe4]*rate_eval.He4_O18__Ne22
       )

    jac[jne22, jf19] = (
       +rho*Y[jhe4]*rate_eval.He4_F19__p_Ne22
       )

    jac[jne22, jne22] = (
       -rate_eval.Ne22__He4_O18
       -rho*Y[jp]*rate_eval.p_Ne22__Na23
       -rho*Y[jp]*rate_eval.p_Ne22__He4_F19
       )

    jac[jne22, jna22] = (
       +rate_eval.Na22__Ne22__weak__wc12
       )

    jac[jne22, jna23] = (
       +rate_eval.Na23__p_Ne22
       )

    jac[jna21, jp] = (
       +rho*Y[jne20]*rate_eval.p_Ne20__Na21
       )

    jac[jna21, jhe4] = (
       +rho*Y[jf17]*rate_eval.He4_F17__Na21
       )

    jac[jna21, jf17] = (
       +rho*Y[jhe4]*rate_eval.He4_F17__Na21
       )

    jac[jna21, jne20] = (
       +rho*Y[jp]*rate_eval.p_Ne20__Na21
       )

    jac[jna21, jna21] = (
       -rate_eval.Na21__Ne21__weak__wc12
       -rate_eval.Na21__p_Ne20
       -rate_eval.Na21__He4_F17
       )

    jac[jna22, jp] = (
       +rho*Y[jne21]*rate_eval.p_Ne21__Na22
       )

    jac[jna22, jhe4] = (
       +rho*Y[jf18]*rate_eval.He4_F18__Na22
       )

    jac[jna22, jf18] = (
       +rho*Y[jhe4]*rate_eval.He4_F18__Na22
       )

    jac[jna22, jne21] = (
       +rho*Y[jp]*rate_eval.p_Ne21__Na22
       )

    jac[jna22, jna22] = (
       -rate_eval.Na22__Ne22__weak__wc12
       -rate_eval.Na22__p_Ne21
       -rate_eval.Na22__He4_F18
       )

    jac[jna23, jp] = (
       -rho*Y[jna23]*rate_eval.p_Na23__He4_Ne20
       +rho*Y[jne22]*rate_eval.p_Ne22__Na23
       )

    jac[jna23, jhe4] = (
       +rho*Y[jf19]*rate_eval.He4_F19__Na23
       +rho*Y[jne20]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jna23, jf19] = (
       +rho*Y[jhe4]*rate_eval.He4_F19__Na23
       )

    jac[jna23, jne20] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne20__p_Na23
       )

    jac[jna23, jne22] = (
       +rho*Y[jp]*rate_eval.p_Ne22__Na23
       )

    jac[jna23, jna23] = (
       -rate_eval.Na23__p_Ne22
       -rate_eval.Na23__He4_F19
       -rho*Y[jp]*rate_eval.p_Na23__He4_Ne20
       )

    return jac
