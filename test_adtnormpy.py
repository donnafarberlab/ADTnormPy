# You must create an editable install:
# python -m pip install -e .
# python -m pytest --cov=adtnormpy --cov-report=html 


import os
os.environ["R_HOME"] = f"{os.environ['CONDA_PREFIX']}/lib/R"
print(os.environ["R_HOME"])

import pytest
import adtnormpy
import pandas as pd
import numpy as np
import rpy2.robjects.packages
np.random.seed(42)

# Testing is currently incomplete. To add:
# Test all of ADTnorm kwargs can be processed without error
# Test on single dataset single marker, multi dataset single marker, and many datasets single marker. Autogenerate gaussians to work with here.
# Write TODO functions for code that we haven't been able to test yet. and test out all kwargs options. 
# Test input as anndata, mudata, or DataFrames 
# Test input errors to make sure they hit assertion errors.
# Add other input checks, e.g. same length between obs and .X
# TODO AssertionError: Ran into R Runtime Error replace RRuntimeError

def assert_order_preserved(data,new_data,obs,include_zeroes = True):      
    assert len(data) == len(new_data)
    assert all(data.index == new_data.index)
    for i in new_data.columns:
        for j in obs['sample'].unique():
            mask = (data[i] != 0) if not include_zeroes else (data[i] != np.nan)
            data1 = data[(obs['sample']==j) & (mask)]
            df1 = new_data[(obs['sample']==j) & (mask)]
            assert all((sorted(df1.loc[data1.sort_values(i).index][i]) == df1.loc[data1.sort_values(i).index][i])), f'Order of {j} donor not preserved {i}'
    return
    
def test_r_import():
    ADTnormR = adtnormpy.ADTnormPy.load_adtnorm_R()
    assert type(ADTnormR) is rpy2.robjects.packages.InstalledSTPackage
    
def test_empty():
    with pytest.raises(TypeError):
        adtnormpy.adtnorm()
        
def test_invalid_dtype():
    with pytest.raises(AssertionError):
        adtnormpy.adtnorm(np.array([1,2,3]))
        
np.random.RandomState(seed=42)
one_data = pd.DataFrame(dict(one_peak=np.concatenate([np.random.normal(loc=1, scale=2.0, size=100),
                                                  np.random.normal(loc=2, scale=2.0, size=100)])))
obs = pd.DataFrame(dict(sample=['1']*100+['2']*100,batch='2'),index=one_data.index)        

def test_just_obs():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(pd.DataFrame(),obs=obs,sample_column='sample',ADT_location=None,save_outpath='tmp/ADTnorm')

def test_no_obs():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(one_data,None,'sample','sample',ADT_location=None,save_outpath='tmp/ADTnorm')

def test_non_counts():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(one_data[one_data.one_peak>0],obs[one_data.one_peak>0],'sample','sample',ADT_location=None,save_outpath='tmp/ADTnorm')
    

int_data = (one_data*300).astype(int)
int_data[int_data.one_peak<30] = 0

def test_two_samples():
    print(int_data)
    print(obs)
    df = adtnormpy.adtnorm(int_data,obs=obs,sample_column='sample',batch_column='sample',ADT_location=None,save_outpath='tmp/ADTnorm')
    assert all(df.columns == ['one_peak'])
    assert_order_preserved(int_data,df,obs,include_zeroes = True)
    
def test_invalid_kwarg_dtype():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(int_data,obs=obs,sample_column='sample',ADT_location=None,save_outpath='tmp/ADTnorm',trimodal_marker =anndata.AnnData())

def test_invalid_kwarg_list_dtype():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(int_data,obs=obs,sample_column='sample',ADT_location=None,save_outpath='tmp/ADTnorm',trimodal_marker =[anndata.AnnData(),anndata.AnnData()])
    
data_scramble = int_data.sample(frac=1,random_state=42)
obs_scramble = obs.loc[data_scramble.index].copy()

def test_scrambled():
    df = adtnormpy.adtnorm(data_scramble,sample_column='sample',obs=obs_scramble,ADT_location=None,save_outpath='tmp/ADTnorm')
    assert all(df.columns == ['one_peak'])
    assert_order_preserved(data_scramble,df,obs_scramble,include_zeroes = True)
    
data = pd.DataFrame(dict(one_peak=np.concatenate([np.random.normal(loc=1, scale=2.0, size=100),
                                                  np.random.normal(loc=2, scale=2.0, size=100)]),
                         one_peak2=np.concatenate([np.random.normal(loc=2, scale=1.0, size=100),
                                                  np.random.normal(loc=6, scale=2.0, size=100)])))
data.loc[data.one_peak < .1,'one_peak'] = 0
data.loc[data.one_peak2 < .1,'one_peak2'] = 0

def test_two_not_raw():    
    df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',save_outpath='tmp/ADTnorm')
    assert all(df.columns == ['one_peak','one_peak2'])
    assert_order_preserved(data,df,obs,include_zeroes = True)
    
def test_two_not_raw_exclude_zeroes():    
    df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',save_outpath='tmp/ADTnorm',exclude_zeroes=True)
    assert all(df.columns == ['one_peak','one_peak2'])
    assert_order_preserved(data,df,obs,include_zeroes = False)
    
def test_view():
    df = adtnormpy.adtnorm(data[data.one_peak>0.5],obs=obs[data.one_peak>0.5],sample_column='sample',ADT_location=None,save_fig=False,save_landmark=False)
    assert_order_preserved(data[data.one_peak>0.5],df,obs[data.one_peak>0.5],include_zeroes = True)

def test_no_save_outpath():
    df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',ADT_location=None,save_outpath=None, save_fig=False,save_landmark=False)
    assert_order_preserved(data,df,obs,include_zeroes = True)
    
def test_no_save_outpath_yes_save_intermediate():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',ADT_location=None,save_outpath=None, save_fig=True,save_landmark=True)
        
# def test_target_landmark_location():
#     df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',ADT_location=None,save_outpath=None,save_fig=False,save_landmark=False,target_landmark_location=(1,5))
#     assert all(df.columns == ['one_peak','one_peak2'])
#     assert_order_preserved(data,df,obs,include_zeroes = True)
    
def test_negative_target_landmark_location():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',ADT_location=None,save_outpath=None,save_fig=False,save_landmark=False,target_landmark_location=(-3,1))
        
def test_single_target_landmark_location():
    with pytest.raises(AssertionError):
        df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',ADT_location=None,save_outpath=None,save_fig=False,save_landmark=False,target_landmark_location=1)

def test_process_one():
    df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',marker_to_process='one_peak',ADT_location=None,save_outpath=None,save_fig=False,save_landmark=False)
    assert all(df.columns == ['one_peak'])
    
def test_process_one_list():
    df = adtnormpy.adtnorm(data,obs=obs,sample_column='sample',marker_to_process=['one_peak'],ADT_location=None,save_outpath=None,save_fig=False,save_landmark=False)
    assert all(df.columns == ['one_peak'])
    assert_order_preserved(data,df,obs,include_zeroes = True)
    
all_modal = pd.DataFrame(dict(one_peak=np.concatenate([np.random.normal(loc=4, scale=.5, size=100),
                                                      np.random.normal(loc=8, scale=2, size=100),
                                                      np.random.normal(loc=9, scale=5, size=100),
                                                      np.random.normal(loc=1, scale=.1, size=100)]),
                             bimodal=np.concatenate([np.random.normal(loc=2, scale=.5, size=40), 
                                                         np.random.normal(loc=5, scale=.5, size=60),
                                                     np.random.normal(loc=3, scale=.5, size=50), 
                                                         np.random.normal(loc=8, scale=.25, size=50),
                                                     np.random.normal(loc=1, scale=.5, size=20), 
                                                         np.random.normal(loc=7, scale=.75, size=80),
                                                     np.random.normal(loc=2, scale=.5, size=60), 
                                                         np.random.normal(loc=5, scale=.5, size=40)]), 
                             trimodal=np.concatenate([np.random.normal(loc=2, scale=.5, size=40), 
                                                         np.random.normal(loc=5, scale=.5, size=40),
                                                         np.random.normal(loc=8, scale=.5, size=20),
                                                     np.random.normal(loc=2, scale=.5, size=30), 
                                                         np.random.normal(loc=5, scale=.5, size=5),
                                                         np.random.normal(loc=9, scale=.75, size=65),                                                    
                                                     np.random.normal(loc=2, scale=.5, size=60), 
                                                         np.random.normal(loc=6, scale=.3, size=10),
                                                         np.random.normal(loc=10, scale=.5, size=30),
                                                    np.random.normal(loc=2, scale=.5, size=33), 
                                                         np.random.normal(loc=5, scale=.5, size=33),
                                                         np.random.normal(loc=8, scale=.5, size=34)])))
all_modal[all_modal<0.1] = 0
all_modal.index = all_modal.index.astype(str)
obs2 = pd.DataFrame(dict(sample=['1']*100+['2']*100+['3']*100+['4']*100,batch='2'),index=all_modal.index)        
all_modal = all_modal.sample(frac=1,random_state=42)
obs2 = obs2.loc[all_modal.index].copy()

def test_all_modal():
    df = adtnormpy.adtnorm(all_modal,obs=obs2,sample_column='sample',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
    assert all(df.columns == ['one_peak','bimodal','trimodal'])
    assert_order_preserved(all_modal,df,obs2,include_zeroes = False)

def test_all_modal_verbose():
    df = adtnormpy.adtnorm(all_modal,obs=obs2,sample_column='sample',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal',verbose=3)
    assert all(df.columns == ['one_peak','bimodal','trimodal'])
    assert_order_preserved(all_modal,df,obs2,include_zeroes = False)
    
def test_all_modal_clean_adt_names():
    df = adtnormpy.adtnorm(all_modal,obs=obs2,sample_column='sample',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal=['bimodal'],trimodal=['trimodal'],clean_adt_name=True)
    assert all(df.columns == ['one','bimodal','trimodal'])
    
import anndata 
adata_in_X = anndata.AnnData(all_modal.astype('float32').values,var=pd.DataFrame(index=all_modal.columns),obs=obs2)

def test_adata_in_X_wrong():
    with pytest.raises(AssertionError):
        data = adata_in_X
        new_data = adtnormpy.adtnorm(data,obs=None,sample_column='sample',ADT_location='protein',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
        
def test_adata_in_X():
    data = adata_in_X
    df = adtnormpy.adtnorm(data,obs=None,sample_column='sample',ADT_location=None,save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
    assert 'ADTnorm' in df.layers.keys()
    assert np.allclose(data.X, df.X)
    obs = data.obs
    data = pd.DataFrame(data.X,index=data.obs_names,columns=data.var_names)
    df = pd.DataFrame(df.layers['ADTnorm'],index=df.obs_names,columns=df.var_names)
    assert_order_preserved(data,df, obs,False)
    
def test_adata_in_X_loc_none():
    data = adata_in_X.copy()
    df = adtnormpy.adtnorm(data,obs=None,sample_column='sample',ADT_location=None,return_location=None,save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
    assert 'ADTnorm' in df.layers.keys()
    assert np.allclose(data.X, df.X)
    obs = data.obs
    data = pd.DataFrame(data.X,index=data.obs_names,columns=data.var_names)
    df = pd.DataFrame(df.X,index=df.obs_names,columns=df.var_names)
    assert_order_preserved(data,df, obs,False)
    
      
adata_in_layers = anndata.AnnData(np.zeros(all_modal.shape)-1,var=pd.DataFrame(index=all_modal.columns),obs=obs2)
adata_in_layers.layers['protein'] = all_modal.values

def test_adata_in_layers():
    data = adata_in_layers
    df = adtnormpy.adtnorm(data,obs=None,sample_column='sample',ADT_location='protein',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
    assert 'ADTnorm' in df.layers.keys()
    assert np.allclose(data.X, df.X)
    assert np.allclose(data.layers['protein'], df.layers['protein'])
    obs = data.obs
    data = pd.DataFrame(data.layers['protein'],index=data.obs_names,columns=data.var_names)
    df = pd.DataFrame(df.layers['ADTnorm'],index=df.obs_names,columns=df.var_names)
    assert_order_preserved(data,df, obs,False)        

adata_in_obsm = anndata.AnnData(np.zeros(all_modal.shape)-1,obs=obs2)
all_modal_prot = all_modal.copy()
all_modal_prot.index = adata_in_obsm.obs_names
adata_in_obsm.obsm['protein'] = all_modal_prot
    
def test_adata_in_obsm():
    data = adata_in_obsm
    df = adtnormpy.adtnorm(data,obs=None,sample_column='sample',ADT_location='protein',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
    assert 'ADTnorm' in df.obsm.keys()
    assert np.allclose(data.X, df.X)
    assert all(data.obsm['protein'] == df.obsm['protein'])
    obs = data.obs
    data = data.obsm['protein']
    df = df.obsm['ADTnorm']
    assert_order_preserved(data,df, obs,False)
    

    
def test_mudata():
    import mudata
    mdata = mudata.MuData(dict(protein=adata_in_X))
    
    data = mdata
    df = adtnormpy.adtnorm(data,obs=None,sample_column='sample',ADT_location='protein',save_outpath='tmp/ADTnorm',exclude_zeroes=True,bimodal='bimodal',trimodal='trimodal')
    assert 'ADTnorm' in df.mod.keys()
    assert np.allclose(data.mod['protein'].X, df.mod['protein'].X)
    obs = data.mod['protein'].obs
    data = pd.DataFrame(data.mod['protein'].X,index=data.obs_names,columns=data.var_names)
    df = pd.DataFrame(df.mod['ADTnorm'].X,index=df.obs_names,columns=df.var_names)
    assert_order_preserved(data,df, obs,False)        

# TODO, add in testing for override_landmark as dict, rds, or csv path, add in tests for landmarks_to_r, landmarks_to_python, save_python_landmarks, and load_python_landmarks