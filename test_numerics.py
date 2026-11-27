"""Small mathematical and protocol checks: python -m unittest -v."""
import unittest
import copy
import numpy as np
import pandas as pd
import analysis as a

class NumericalChecks(unittest.TestCase):
    def setUp(self):
        self.cfg=a.configuration()
        self.x=np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.],[2.,1.]])
        self.y=np.array([1.,2.,3.,4.,5.]); self.e=np.zeros(5)
    def test_idw_exact_and_lapse(self):
        p=a.predict('IDW',self.x,self.y,self.x,self.e,self.e,self.cfg)
        np.testing.assert_allclose(p,self.y)
        p=a.predict('IDW',self.x,self.y,self.x,self.e,np.ones(5)*100,self.cfg)
        np.testing.assert_allclose(p,self.y-.65)
    def test_kriging_exact_unbiased(self):
        w=a.kriging_weights(self.x,self.x,0,1,2)
        np.testing.assert_allclose(self.y@w,self.y,atol=1e-10)
        np.testing.assert_allclose(w.sum(axis=0),1)
    def test_constant_temperature_still_uses_idw_elevation(self):
        p=a.predict('IDW',self.x,np.ones(5),self.x,self.e,np.ones(5)*100,self.cfg)
        np.testing.assert_allclose(p,.35)
    def test_duplicate_rejected(self):
        with self.assertRaises(ValueError): a.validate_xy(np.zeros((3,2)),np.ones(3),self.x)
    def test_metrics(self):
        m=a.metrics([1,2,3],[2,2,2])
        self.assertAlmostEqual(m['MAE'],2/3); self.assertAlmostEqual(m['MSE'],2/3); self.assertAlmostEqual(m['NSE'],0)
        self.assertTrue(np.isnan(a.metrics([1,1],[2,2])['NSE']))
        with self.assertRaises(ValueError): a.metrics([1,2],[1,np.nan])
        with self.assertRaises(ValueError): a.metrics([1,2],[[1],[2]])
    def test_normal_ties_inverse(self):
        y=np.array([1.,1.,2.,5.]); z=a.normal_scores(y)
        self.assertEqual(z[0],z[1]); np.testing.assert_allclose(a.inverse_scores(z,y),y)
    def test_topsis_negative_constant(self):
        df=pd.DataFrame(dict(method=['a','b','c'],MAE=[1,2,3],NSE=[-1,-2,-3],constant=[2,2,2]))
        r,w=a.topsis(df,dict(MAE=-1,NSE=1,constant=1))
        self.assertEqual(w['constant'],0); self.assertEqual(r['rank'].tolist(),[1,2,3])
        df['MAE']=1; df['NSE']=0
        r,_=a.topsis(df,dict(MAE=-1,NSE=1)); np.testing.assert_allclose(r.closeness,.5)
    def test_sampling_bounded_reproducible(self):
        h,g=a.empirical_variogram(self.x,a.normal_scores(self.y),10)
        s,_=a.sample_parameters(h,g,self.cfg,np.random.default_rng(4))
        t,_=a.sample_parameters(h,g,self.cfg,np.random.default_rng(4))
        np.testing.assert_array_equal(s,t)
        bounds=np.array(self.cfg['parameter_bounds']); self.assertTrue(((s>=bounds[:,0])&(s<=bounds[:,1])).all())
    def test_local_parallel_and_batch_reproducible(self):
        q=np.array([[.2,.1],[.7,.6]])
        p=a.predict('EBK-inspired local kriging',self.x,self.y,q,self.e,np.zeros(2),self.cfg,4)
        cfg=copy.deepcopy(self.cfg); cfg['workers']=2
        r=a.predict('EBK-inspired local kriging',self.x,self.y,q,self.e,np.zeros(2),cfg,4)
        np.testing.assert_array_equal(p,r)
    def test_loo_withheld_value_excluded(self):
        df=pd.DataFrame(dict(longitude=[49,49.25,49,49.25,49.5],latitude=[37,37,37.25,37.25,37.25],temperature=self.y,elevation=self.e))
        pred,_=a.cross_validate({'annual':df},self.cfg)
        changed=df.copy(); changed.loc[0,'temperature']=999
        again,_=a.cross_validate({'annual':changed},self.cfg)
        np.testing.assert_array_equal(pred[pred.location_id==0].prediction,again[again.location_id==0].prediction)

if __name__=='__main__': unittest.main()
