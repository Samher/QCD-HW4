import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

def dGammaChiral(m12sq):
    return 1 / (m_etaprime**4) * 1/4 * (m12sq - 2*m_pi**2)**2


def int1(m12sq, matel):
    E2 = 1 / 2 * np.sqrt(m12sq)
    E3 = (m_etaprime**2 - m12sq - m_eta**2) / (2*np.sqrt(m12sq))
    m23sqmin = (E2 + E3)**2 - (np.sqrt(E2**2 - m_pi**2) + np.sqrt(E3**2 - m_pi**2))
    m23sqmax = (E2 + E3)**2 - (np.sqrt(E2**2 - m_pi**2) - np.sqrt(E3**2 - m_pi**2))
    return (m23sqmax - m23sqmin) * matel(m12sq)


def int2(matel):
    result = spi.quad(int1, (2 * m_pi)**2, (m_etaprime - m_eta)**2, args=(matel))
    return result[0]


m_etaprime = 957.78  # MeV
m_eta = 547.862      # MeV
m_pi = 139.58039     # MeV

GammaChiral = int2(dGammaChiral)
GammaNaive = int2(lambda _: 1)

print(GammaChiral/GammaNaive)
