#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sklearn.datasets import load_boston
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.impute import SimpleImputer
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   #Data visualisation libraries 
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
get_ipython().run_line_magic('matplotlib', 'inline')


# In[4]:


# Loading Data
boston_data = load_boston()
print(boston_data.data.shape)


# In[5]:


boston = pd.DataFrame(boston_data.data,columns=boston_data.feature_names)
boston['target'] = pd.Series(boston_data.target)
boston.head()


# In[6]:


# In original dataset there is no missings
print('Dataframe Null values:',boston.isnull().values.ravel().sum())


# In[7]:


# Number of unique categories
print ('Number of unique classes:\n',boston.nunique())


# In[8]:


sns.pairplot(boston)


# In[9]:


# Distribution plot of target variable which is MEDV
sns.distplot(boston['target'])


# # Step 1:
# 
# Using Sklearn get the Boston Housing dataset.
# Fit a linear regressor to the data as a baeline.  There is no need to do Cross-Validation.  We are exploring the change in results
# What is the loss and what are the goodness of fit parameters?  This will be our baseline for comparison
# 

# In[10]:


# dividing dataset to two X as explanatory variables and Y as target variable
X=boston.iloc[:,0:13]
Y=boston['target']


# In[11]:


# Fitting a line
model = LinearRegression().fit(X, Y)
baselineR2 = model.score(X, Y) 
print('Baseline R2:',baselineR2)
# the coefficients of different variables
model.coef_


# In[12]:


# graph to compare actual values and predicted values of target
baseline_predict=model.predict(X)
plt.scatter(Y, baseline_predict)


# In[13]:


# Calculating Scores for the baseline regression
# mean square error
base_MSE=mean_squared_error(Y, baseline_predict)
print('Baseline MSE is:',base_MSE)
# mean_absolute_error
base_MAE=mean_absolute_error(Y,baseline_predict)
print('Baseline MAE is:',base_MAE)
# r2_score
base_R2=r2_score(Y, baseline_predict)
print('Baseline R2 score is:',base_R2)


# # Step 2: (repeated)
# For select between 1, 5 10, 20, 33, and 50% of your data on a single column (Completely at random), replace the present value with a NAN and then perform an imputation of that value.   
# 
# In. each case perform a fit with the imputed data and compare the loss and goodness of fit to your baseline.
# 

# In[14]:


def step2_results(dataset,target,percent,imp_var,strategy):
    # Set percentage of values randomly to NA
    if percent != 0:
        np.random.seed(42)
        indexer = np.sort(np.random.permutation(len(dataset))[len(dataset)-(int(len(dataset)*percent)):])
        dataset_imp = dataset.copy()
        # Chose LSTAT varialbe to randomly set to NaN. Chose LSTAT becuase of the 
        # large t-statistic in the original analysis
        dataset_imp[imp_var][indexer] = np.nan
        
        imp_model = SimpleImputer(missing_values=np.nan, strategy=strategy)
        imp_model.fit(dataset_imp)
        
        imp_model = pd.DataFrame(imp_model.transform(dataset_imp),columns=dataset.columns, index=dataset.index)
    else:
        imp_model = dataset.copy()
    
    newX = sm.add_constant(imp_model)
    regOLS_imp = sm.OLS(target,newX).fit()
    target_pred = model.predict(imp_model)
    mse = mean_squared_error(target, target_pred)
    mae=mean_absolute_error(target, target_pred)
    #r2=r2_score(target, target_pred)
    return((percent,regOLS_imp._results.rsquared,regOLS_imp._results.rsquared_adj,
            regOLS_imp._results.bic,mse,mae,imp_model))


# In[15]:


results_boston = pd.DataFrame([])
per_num = [0.0,0.01,0.05,0.10,0.20,0.33,0.50]
imp_var = ['LSTAT'] #Choose from CRIM,ZN,INDUS,CHAS,NOX,RM,AGE,DIS,RAD,TAX,PTRATIO,B,LSTAT
strategy = 'mean' #Choose from 'mean', 'median', 'constant'
#fill_value = 0 # If you choose 'constant', need to select value to fill NaNs


# In[16]:


for i in range(len(per_num)):
    results_boston[i] = np.array(step2_results(X,Y,per_num[i],imp_var[0],strategy), dtype=object)
results_boston


# In[17]:


results_boston = results_boston.T
results_boston.columns=['Imputed_Percent','RSquared','AdjRSquared','BIC','MSE','MAE','Model']

print(results_boston.loc[:,results_boston.columns !='Model']) # Prints results of analysis for each percent of Imputed variables


# In[18]:


# Plot of R2 values analysis results vs different Percentage of Imputation
results_boston.plot('Imputed_Percent','RSquared')
plt.ylabel('RSquared Value')
plt.title('RSquared Value vs Imputed Percentage')

# Plot of MSE values analysis results vs different Percentage of Imputation
results_boston.plot('Imputed_Percent','MSE')
plt.ylabel('Mean Square Error')
plt.title('Mean Square Error vs Imputed Percentage')

# Plot of MAE values analysis results vs different Percentage of Imputation
results_boston.plot('Imputed_Percent','MAE')
plt.ylabel('Mean Absolute Error')
plt.title('Mean Absolute Error vs Imputed Percentage')


# In[19]:


imputed_dataframe = pd.DataFrame([])
plot_result = pd.DataFrame([])
for k in range(0,len(results_boston)):
    imputed_dataframe = results_boston.Model[k]
    imputed_dataframe['Imputed_Percent'] = per_num[k]
    plot_result = plot_result.append(imputed_dataframe, ignore_index=True)

#  Boxplot of Imputed Variable vs different Percentage of Imputation
sns.set(style="whitegrid")
plt.figure(figsize=(16,10))
sns.set(font_scale=1.5)
ax = sns.boxplot(x="Imputed_Percent", y="LSTAT", data=plot_result)
ax.set_title('Box Plot - LSTAT Variable vs Percent Imputed')


summaryResult = pd.DataFrame([])
for k in range(0,len(results_boston)):
    tempResults = results_boston.Model[k]['LSTAT'].describe()
    tempResults['Imputed_Percent']=per_num[k]
    summaryResult = summaryResult.append(tempResults, ignore_index=True)

print(summaryResult)


# # Step 3: 
# Take 2 different columns and create data “Missing at Random” when controlled for a third variable (i.e if Variable Z is > 30, than Variables X, Y are randomly missing).  Make runs with 10%, 20% and 30% missing data imputed via your best guess.  Repeat your fit and comparisons to the baseline.

# In[20]:


# Create a function that perform imputation and fit the model and return the results
def step3_results(dataset,target,percent,ran_seed,cond_var,cond_var_val,imp_var,strategy,fill_value=None):
    datasetDesc = dataset.describe()
    
    if percent != 0:
        Prob3SubIndex = dataset[dataset[cond_var] > datasetDesc.loc[cond_var_val][cond_var]].index
        lenIndex = len(dataset[dataset[cond_var] > datasetDesc.loc[cond_var_val][cond_var]])
        from sklearn.impute import SimpleImputer
        for i in range(len(imp_var)):   
            np.random.seed(ran_seed[i])
            indexer = np.sort(np.random.permutation(Prob3SubIndex)[lenIndex-(int(lenIndex*percent)):])
            step3[imp_var[i]][indexer] = np.nan
            
            imp_modelP3 = SimpleImputer(missing_values=np.nan, strategy=strategy,fill_value=fill_value)
            imp_modelP3.fit(dataset)
                
            imp_modelP3 = pd.DataFrame(imp_modelP3.transform(dataset),columns=dataset.columns, index=dataset.index)
    else:
        imp_modelP3 = dataset.copy()

    X = sm.add_constant(imp_modelP3)
    regOLS_impP3 = sm.OLS(target,X).fit()
    regOLS_impP3.summary()
    regOLS_impP3._results.rsquared
    target_pred = model.predict(imp_modelP3)
    mse = mean_squared_error(target, target_pred) # Not sure about MSE values
    mae=mean_absolute_error(target, target_pred)
    #r2=r2_score(target, target_pred)
    return((percent,regOLS_impP3._results.rsquared,regOLS_impP3._results.rsquared_adj,
            regOLS_impP3._results.bic,mse,mae))


# In[21]:


per_numP3 = [0.0,0.10,0.20, 0.30]  # Percentage of imputed values in dataset
ran_seed = [42,39] #Random seed numbers for permutation
cond_var = 'B'   #Choose from CRIM,ZN,INDUS,CHAS,NOX,RM,AGE,DIS,RAD,TAX,PTRATIO,B,LSTAT
cond_var_val ='25%'   # Choose from 25%, 50%, 75%, mean,min
imp_var = ['LSTAT','PTRATIO'] #Choose from CRIM,ZN,INDUS,CHAS,NOX,RM,AGE,DIS,RAD,TAX,PTRATIO,B,LSTAT
strategy = 'constant' #Choose from 'mean', 'median', 'constant'
fill_value = 0 # If you choose 'constant', need to select value to fill NaNs


# In[22]:


step3 = X.copy() # Making copy of original dataset
#Prob3.describe()
results_bostonP3 = pd.DataFrame([]) # Initalizing results dataframe
for j in range(len(per_numP3)):
    results_bostonP3[j] = np.array(step3_results(step3,Y,per_numP3[j],ran_seed,cond_var,cond_var_val,imp_var,strategy,fill_value))

results_bostonP3 = results_bostonP3.T
results_bostonP3.columns=['Imputed_Percent','RSquared','AdjRSquared','BIC','MSE', 'MAE']

print(results_bostonP3) # Prints results of analysis for each percent of Imputed variables


# In[23]:


results_bostonP3.plot('Imputed_Percent','RSquared')
plt.ylabel('RSquared Value')
plt.title('RSquared Value vs Imputed Percentage')

# Plot of MSE values analysis results vs different Percentage of Imputation
results_bostonP3.plot('Imputed_Percent','MSE')
plt.ylabel('Mean Square Error')
plt.title('Mean Square Error vs Imputed Percentage')

# Plot of MAE values analysis results vs different Percentage of Imputation
results_bostonP3.plot('Imputed_Percent','MAE')
plt.ylabel('Mean Absolute Error')
plt.title('Mean Absolute Error vs Imputed Percentage')


# In[ ]:





# In[ ]:





# # Step 4:  
# Create a Missing Not at Random pattern in which 25% of the data is missing for a single column.    Impute your data, fit the results and compare to a baseline.

# In[24]:


def step4_results(dataset,target,percent,imp_var,strategy,fill_value=None):
    # Set 25% percent of values to NA, selecting every 4th entry
    dataset_imp = dataset.copy()
    indexer = dataset_imp.iloc[::4, :].index
    dataset_imp[imp_var][indexer] = np.nan
    
    
    imp_model = SimpleImputer(missing_values=np.nan, strategy=strategy,fill_value=fill_value)
    imp_model.fit(dataset_imp)
    
    imp_model = pd.DataFrame(imp_model.transform(dataset_imp),columns=dataset.columns, index=dataset.index)
    
    X = sm.add_constant(imp_model)
    regOLS_imp = sm.OLS(target,X).fit()
    #regOLS_imp.summary()
    #regOLS_imp._results.rsquared
    target_pred = model.predict(imp_model)
    mse = mean_squared_error(target, target_pred) 
    mae=mean_absolute_error(target, target_pred)
    #r2=r2_score(target, target_pred)
    return((percent,regOLS_imp._results.rsquared,regOLS_imp._results.rsquared_adj,
            regOLS_imp._results.bic,mse,mae))


# In[25]:


# our baseline was done with "LSTAT" column and in order to compare results, we do it the same column
imp_var = ['LSTAT']
strategy = 'mean' 
fill_value = 0 # If you choose 'constant', need to select value to fill NaNs
per_numP4 = [0.25]  # 25% of imputation


# In[26]:


step4 = X.copy() # Making copy of original dataset
#Prob4.describe()
results_bostonP4 = pd.DataFrame([]) # Initalizing results dataframe
for j in range(len(imp_var)):
    results_bostonP4[j] = np.array(step4_results(step4,Y,per_numP4[j],imp_var[j],strategy,fill_value))

results_bostonP4 = results_bostonP4.T
results_bostonP4.columns=['PercentImputed','RSquared','AdjRSquared','BIC','MSE','MAE']

print(results_bostonP4) # Prints results of analysis for each percent of Imputed variables


# In[27]:


# Now we just need to compare above results with the step1 results (baseline)


# In[28]:


base_df = pd.DataFrame([])
base_df = base_df.append({'PercentImputed': 0, 'RSquared': baselineR2, 'MSE': base_MSE, 'MAE': base_MAE}, ignore_index=True)

baseline_comparison_df = pd.DataFrame(columns=results_bostonP4.columns)
baseline_comparison_df = baseline_comparison_df.append(base_df)
baseline_comparison_df = baseline_comparison_df.append(results_bostonP4)
print(baseline_comparison_df)

# Now we just need to compare above results with the step1 results (baseline)
# Plot of R2 values analysis results vs different Percentage of Imputation

baseline_comparison_df.plot('PercentImputed', 'RSquared')
plt.ylabel('RSquared Value')
plt.title('RSquared Value vs Imputed Percentage')

# Plot of MSE values analysis results vs different Percentage of Imputation
baseline_comparison_df.plot('PercentImputed', 'MSE')
plt.ylabel('Mean Square Error')
plt.title('Mean Square Error vs Imputed Percentage')

# Plot of MAE values analysis results vs different Percentage of Imputation
baseline_comparison_df.plot('PercentImputed', 'MAE')
plt.ylabel('Mean Absolute Error')
plt.title('Mean Absolute Error vs Imputed Percentage')


# # Step 5 (Extra Credit) (10 points):
# Using the MCMC method, and your data from step 4, What is the difference in performance between imputation via ‘guess’ (mean/median, etc) and MCMC. 

# In[38]:



dataset = X.copy()
imputed_variable_name = 'LSTAT'
imputed_dataset = dataset.copy()
indexer = imputed_dataset.iloc[::4, :].index
imputed_dataset[imputed_variable_name][indexer] = np.nan
sns.set(color_codes=True)
sns.distplot(imputed_dataset['LSTAT'], hist=False, rug=True);


# In[ ]:




