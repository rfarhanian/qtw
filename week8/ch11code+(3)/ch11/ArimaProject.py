# matplotlib inline
import datetime
import datetime
import numpy as np
import pandas as pd
from statsmodels.tsa.stattools import acf, pacf
from statsmodels.tsa.arima_model import ARIMA
from pandas_datareader import data as web
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.tsa.stattools import adfuller

pd.core.common.is_list_like = pd.api.types.is_list_like
from pandas_datareader import data as web

start = datetime.datetime(2016, 1, 1)
end = datetime.datetime(2019, 2, 1)
# price = web.DataReader('BABA', 'yahoo', start, end)
# price.index = pd.to_datetime(price.index)
#
#
# # data = read_csv('AirPassengers.csv',sep=',')
# data = web.DataReader('BABA', 'yahoo', start, end)
# stock_price = data['Close']
# data.index = pd.to_datetime(data.index)
# print(stock_price.head())
# fig = plt.figure()
#
# stock_price.plot()
# plt.title("AliBaba Stocks 2016-2019")
# plt.xlabel("Time")
# plt.ylabel('Price')
# plt.show()
#
# data.describe
#

#
#
def test_stationarity(timeseries):
    # Determing rolling statistics
    rolmean = timeseries.rolling(12).mean()
    rolstd = timeseries.rolling(12).std()

    # Plot rolling statistics:
    orig = plt.plot(timeseries, color='blue', label='Original')
    mean = plt.plot(rolmean, color='red', label='Rolling Mean')
    std = plt.plot(rolstd, color='black', label='Rolling Std')
    plt.legend(loc='best')
    plt.title('Rolling Mean & Standard Deviation')
    plt.show(block=False)

    # Perform Dickey-Fuller test:
    print('Results of Dickey-Fuller Test:')
    dftest = adfuller(timeseries, autolag='AIC')
    dfoutput = pd.Series(dftest[0:4], index=['Test Statistic', 'p-value', '#Lags Used', 'Number of Observations Used'])
    for key, value in dftest[4].items():
        dfoutput['Critical Value (%s)' % key] = value
    print(dfoutput)
#
#
# test_stationarity(stock_price)
#
# moving_avg = stock_price.rolling(12).mean()
# no_trend = stock_price - moving_avg
# plt.plot(no_trend)
# plt.title("Stock Data with Rolling Mean Removed")
# test_stationarity(no_trend.dropna())
#
# log_stock_price = np.log(stock_price)
# log_moving_avg = log_stock_price.rolling(12).mean()
# log_no_trend = log_stock_price - log_moving_avg
# test_stationarity(log_no_trend.dropna())
#
# from statsmodels.tsa.seasonal import seasonal_decompose
#
# decomposition = seasonal_decompose(stock_price)
#
# trend = decomposition.trend
# seasonal = decomposition.seasonal
# residual = decomposition.resid
# plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
# plt.subplot(411)
# plt.plot(stock_price, label='Original')
# plt.legend(loc='best')
# plt.subplot(412)
# plt.plot(trend, label='Trend')
# plt.legend(loc='best')
# plt.subplot(413)
# plt.plot(seasonal, label='Seasonality')
# plt.legend(loc='best')
# plt.subplot(414)
# plt.plot(residual, label='Residuals')
# plt.legend(loc='best')
# plt.tight_layout()
#
# diff_stock_price = log_stock_price - log_stock_price.shift()
# test_stationarity(diff_stock_price.dropna())
#
# acf_stock_price = acf(diff_stock_price.dropna(), nlags=20)
# # Plot ACF:
# plt.figure(figsize=(20, 4))
#
# plt.subplot(121)
# plt.plot(acf_stock_price)
# plt.xticks(np.arange(21))
# plt.axhline(y=0, linestyle='--', color='gray')
# plt.axhline(y=-1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
# plt.axhline(y=1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
# for i in range(1, 20):
#     plt.axvline(x=i, linestyle=':', color='gray')
# plt.title('Autocorrelation Function')
#
# pacf_air_plot = pacf(diff_stock_price.dropna(), nlags=20)
# # Plot ACF:
# plt.figure(figsize=(20, 4))
#
# plt.subplot(121)
# plt.plot(pacf_air_plot)
# plt.xticks(np.arange(21))
# plt.axhline(y=0, linestyle='--', color='gray')
# plt.axhline(y=-1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
# plt.axhline(y=1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
# for i in range(1, 20):
#     plt.axvline(x=i, linestyle=':', color='gray')
# plt.title('Partial Autocorrelation Function')
#
#
#
# model = ARIMA(diff_stock_price.dropna(), order=(1, 0, 1))
# model_fit = model.fit(disp=0)
# print(model_fit.summary())
#
# # plot residual errors
#
# residuals = pd.DataFrame(model_fit.resid)
# residuals.plot()
# plt.show()
# residuals.plot(kind='kde')
# plt.show()
# print(residuals.describe())
#
# stock_price[-49:].plot()
#
# stock_price = stock_price.astype(float)
# loss_best = 1E16
# best_ints = [-1, -1, -1]
# for p in range(4):
#     for d in range(3):
#         for q in range(3):
#             model = ARIMA(diff_stock_price.dropna(), order=(p, d, q))
#             try:
#                 results_ARIMA = model.fit(disp=-1)
#             except ValueError:
#                 pass
#             except:
#                 pass
#             plt.plot(diff_stock_price)
#             plt.plot(results_ARIMA.fittedvalues, color='red')
#             x = pd.DataFrame(results_ARIMA.fittedvalues)
#             x = x.join(diff_stock_price)
#             x['out'] = (x.iloc[:, 0] - x.iloc[:, 1]) ** 2
#             loss = np.sqrt(x['out'].sum())
#             plt.title('RSS: %.4f' % loss)
#             if loss < loss_best:
#                 print(loss)
#                 loss_best = loss
#                 best_ints = [p, d, q]
#             plt.show()
#             print(p, d, q)
#
# stock_price = stock_price.astype(float)
# loss_best = 1E16
# best_ints = [-1, -1, -1]
# for p in range(4):
#     for d in range(3):
#         for q in range(3):
#             model = ARIMA(stock_price.dropna(), order=(p, d, q))
#             try:
#                 results_ARIMA = model.fit(disp=-1)
#             except ValueError:
#                 pass
#             except:
#                 pass
#             plt.plot(stock_price)
#             plt.plot(results_ARIMA.fittedvalues, color='red')
#             x = pd.DataFrame(results_ARIMA.fittedvalues)
#             x = x.join(stock_price)
#             x['out'] = (x.iloc[:, 0] - x.iloc[:, 1]) ** 2
#             loss = np.sqrt(x['out'].sum())
#             plt.title('RSS: %.4f' % loss)
#             if loss < loss_best:
#                 print(loss)
#                 loss_best = loss
#                 best_ints = [p, d, q]
#             plt.show()
#             print(p, d, q)
#
# print(loss_best)
# print(best_ints)
#
#
#
# pd.core.common.is_list_like = pd.api.types.is_list_like


price = web.DataReader('AAPL', 'yahoo', start, end)['Close']
price.index = pd.to_datetime(price.index)

price.plot()

test_stationarity(price)

diff = price - price.shift()
diff.dropna(inplace=True)
diff.plot()
test_stationarity(diff)

acf_stock_price = acf(diff, nlags=20)
# Plot ACF:
plt.figure(figsize=(20, 4))

plt.subplot(121)
plt.plot(acf_stock_price)
plt.xticks(np.arange(21))
plt.axhline(y=0, linestyle='--', color='gray')
plt.axhline(y=-1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
plt.axhline(y=1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
for i in range(1, 20):
    plt.axvline(x=i, linestyle=':', color='gray')
plt.title('Autocorrelation Function')

pacf_air_plot = pacf(diff, nlags=20)
# Plot ACF:
plt.figure(figsize=(20, 4))

plt.subplot(121)
plt.plot(pacf_air_plot)
plt.xticks(np.arange(21))
plt.axhline(y=0, linestyle='--', color='gray')
plt.axhline(y=-1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
plt.axhline(y=1.96 / np.sqrt(len(acf_stock_price)), linestyle='--', color='gray')
for i in range(1, 20):
    plt.axvline(x=i, linestyle=':', color='gray')
plt.title('Partial Autocorrelation Function')

model = ARIMA(diff, order=(1, 0, 0))
model_fit = model.fit(disp=0)
print(model_fit.summary())
# plot residual errors
residuals = pd.DataFrame(model_fit.resid)
residuals.plot()
plt.show()
residuals.plot(kind='kde')
plt.show()
print(residuals.describe())

loss_best = 1E16
best_ints = [-1, -1, -1]
for p in range(4):
    for d in range(3):
        for q in range(3):
            model = ARIMA(diff, order=(p, d, q))
            try:
                results_ARIMA = model.fit(disp=-1)
            except ValueError:
                pass
            except:
                pass
            plt.plot(diff)
            plt.plot(results_ARIMA.fittedvalues, color='red')
            x = pd.DataFrame(results_ARIMA.fittedvalues)
            x = x.join(diff)
            x['out'] = (x.iloc[:, 0] - x.iloc[:, 1]) ** 2
            loss = np.sqrt(x['out'].sum())
            plt.title('RSS: %.4f' % loss)
            if loss < loss_best:
                print(loss)
                loss_best = loss
                best_ints = [p, d, q]
            plt.show()
            print(p, d, q)

loss_best
