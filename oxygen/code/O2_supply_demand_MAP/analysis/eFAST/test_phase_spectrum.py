"""Independent Fourier identities and rejection of invalid design identities."""
import unittest
import numpy as np
import phase_spectrum as diag


class SpectrumChecks(unittest.TestCase):
    def test_known_low_and_focal_sine_power(self):
        n,m=4097,4;s=2*np.pi*np.arange(n)/n;omega=(n-1)//(2*m)
        y=3*np.sin(s)+2*np.cos(omega*s)
        power,v,low,s1,st,_,_=diag.spectrum(y,m)
        self.assertAlmostEqual(float(v),6.5,places=12)
        self.assertAlmostEqual(float(low),4.5,places=12)
        self.assertAlmostEqual(float(s1),2/6.5,places=12)
        self.assertAlmostEqual(float(st),2/6.5,places=12)
        np.testing.assert_allclose(v,np.var(y),rtol=1e-12)

    def test_vectorized_fft_agrees_with_independent_full_fft(self):
        n=8193;m=4;rng=np.random.default_rng(123)
        y=rng.normal(size=(n,3,2))
        _,v,_,s1,st,omega,cut=diag.spectrum(y,m)
        for j in range(3):
            for k in range(2):
                p=(np.abs(np.fft.fft(y[:,j,k])[1:(n+1)//2])/n)**2
                self.assertAlmostEqual(float(v[j,k]),2*sum(p),places=12)
                self.assertAlmostEqual(float(s1[j,k]),2*sum(p[np.arange(1,m+1)*omega-1])/(2*sum(p)),places=12)
                self.assertAlmostEqual(float(st[j,k]),1-sum(p[:cut])/sum(p),places=12)

    def test_phase_identity_and_wrong_M_rejected(self):
        n,m,p=4097,4,10;phi=4.3
        # Algebraic unit fixture only; not a new scientific sample or model run.
        w=diag.frequencies(n,m,p);s=2*np.pi*np.arange(n)/n
        u=.5+np.arcsin(np.sin(s[:,None]*w+phi))/np.pi
        actual,error,_=diag.infer_phase(u,n,m,p)
        self.assertAlmostEqual(actual,phi,places=10)
        self.assertLess(error,1e-8)
        with self.assertRaises(AssertionError):diag.infer_phase(u,n,8,p)
        with self.assertRaises(AssertionError):diag.infer_phase(u[::-1],n,m,p)

    def test_undefined_or_incomplete_spectrum_rejected(self):
        with self.assertRaises(ValueError):diag.spectrum(np.ones(129),4)
        with self.assertRaises(AssertionError):diag.spectrum(np.full(129,np.nan),4)
        with self.assertRaises(AssertionError):diag.spectrum(np.ones(128),4)


if __name__=='__main__':unittest.main()
