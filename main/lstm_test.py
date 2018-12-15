import pandas as pd
import numpy as np
import matplotlib.dates as mdate
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasRegressor
from keras.layers.normalization import BatchNormalization
from keras.layers import Dropout
#import datetime
from datetime import datetime
from dateutil.relativedelta import relativedelta
import keras
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, DateFormatter
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten, LSTM, TimeDistributed, RepeatVector
from keras.layers.normalization import BatchNormalization
from keras.optimizers import Adam
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping, ModelCheckpoint
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
import csv
import os


def txt2csv(txt_file, csv_file):    
    in_txt = csv.reader(open(txt_file, 'rt'), delimiter = '\t')
    out_csv = csv.writer(open(csv_file, 'w', newline=''))
    out_csv.writerows(in_txt)

def readData(file):
    data = pd.read_csv(file)
    return data
      
def chopDataSetByMonth(data, date):
    for index, row in data.iterrows():
        if date in row['Data']:
            return data[:index], data[index:]
        
def splitData(X,Y, date, rate):
    X_train = X[:int(X.shape[0]*rate)]
    Y_train = Y[:int(Y.shape[0]*rate)]
    train_date = date[:int(X.shape[0]*rate)] 
    X_val = X[int(X.shape[0]*rate):]
    Y_val = Y[int(Y.shape[0]*rate):]
    test_date = date[int(X.shape[0]*rate):]
  
    return X_train, Y_train, X_val, Y_val, train_date, test_date
  
#using max-min normalization  
def MaxMinNormalization(train):
    features = (train.columns.values).tolist()
    features.remove('Data')
    features.remove('exukus')

    #normalize data
    train[features] = train[features].apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)), axis=1)  
    return train
       
# using zscore normalization 
def ZScoreNormalization(train):
    features = (train.columns.values).tolist()
    features.remove('Data')
    features.remove('exukus')  

    train[features] = train[features].apply(lambda x: (x - np.mean(x)) / x.std())   
    return train
 
def buildSimpleRegressionTrainData(train):
    X_train, Y_train, date = [], [], []
    date = train['Data'].values
    Y_train = train['exukus'].values
    X_train = train.drop(['Data', 'exukus'], axis = 1).values
    
    return X_train, Y_train, date

def buildSimpleRegressionTestData(data, sliding_win_size):
    y = data['exukus'].values
    x = data.drop(['Data', 'exukus'], axis = 1).values
    X_test, Y_test, dates = [], [], []
    
    for i in range(data.shape[0]-sliding_win_size + 1):
        X_test.append(x[i:i+sliding_win_size])
        Y_test.append(y[i:i+sliding_win_size])
        dates.append(np.array(data.iloc[i:i+sliding_win_size]["Data"]))
    
    return np.array(X_test), np.array(Y_test), np.array(dates)

def buildTrainData(train, past, future):
    X_train, Y_train, date = [], [], []
    target_column_index = train.columns.get_loc("exukus") -1 
  
    for i in range(train.shape[0]-future-past+1):
        X_train.append(np.array(train.drop(['Data'], axis=1).iloc[i:i+past]))   
        Y_train.append(np.array(train.iloc[i+past:i+past+future]["exukus"]))
        date.append(np.array(train.iloc[i:i+past+future]["Data"]))
            
    return target_column_index, np.array(X_train), np.array(Y_train), np.array(date)

def generateFutureDate(start_date, months_count):
    future_dates = []
    begin = datetime.strptime(start_date, '%Y-%m-%d').date()
    
    for i in range(months_count):
        date = begin + relativedelta(months=+i)
        future_dates.append(date.strftime('%Y-%m-%d'))
    return np.array(future_dates)
        
def buildOnePredictedData(data, future, start_date):
    x = []
    y = []
    date = []
    for index, row in data.iterrows():
        if start_date in row['Data']:
            if (index + future > data.shape[0]):
                print("out of bound")
                return None, None, None
            else:
                x.append(np.array(data.drop(['Data'], axis=1).iloc[index:index+future]))
                if (data.shape[0] - (index + future + future) >= 0):
                    label = np.array(data.iloc[index+future:index+future+future]["exukus"])
                    y.append(label)
                    date.append(np.array(data.iloc[index:index+future + future]["Data"]))
                    return np.array(x), np.array(y), np.array(date)
                else:
                    label = np.array(data.iloc[index+future:]["exukus"])
                    nan_array = np.empty(future - len(label))
                    nan_array.fill(np.nan)
                    y_label = np.concatenate((label, nan_array))
                    date_have = data.iloc[index:]["Data"].values
                    start_date = data.iloc[data.shape[0] - 1]["Data"]
                    future_dates = generateFutureDate(start_date, future - len(label))
                    full_dates = np.concatenate((date_have, future_dates))
                    y.append(y_label)
                    date.append(full_dates)
                    return np.array(x), np.array(y), np.array(date)
                
    print("No such date")
    return None, None, None

def buildSimpleRegressionModel(dim):
    model = Sequential()
    model.add(Dense(30, input_dim=dim, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(1))
    #opt = keras.optimizers.SGD(lr= 0.001)
    opt = keras.optimizers.Adam(lr=0.001)
    model.compile(optimizer=opt, loss='mse', metrics=['mae'])
    model.summary()
    return model

def buildManyToManyModel(shape):
    model = Sequential()
    model.add(LSTM(35, input_length=shape[1], input_dim=shape[2], return_sequences=True))     
    model.add(LSTM(25,return_sequences=True))
    model.add(TimeDistributed(Dense(1)))
    opt = keras.optimizers.RMSprop(lr=0.0001)
    #opt = keras.optimizers.Adam(lr=0.0001)
    model.compile(loss="mse", optimizer=opt, metrics=["mse"])
    model.summary()
    return model

#Many to many model. The length of the input is same as the output
def buildManyToManyModel3(shape, learning_rate=0.0001):
    model = Sequential()
    model.add(LSTM(30, input_length=shape[1], input_dim=shape[2], return_sequences=True))     
    model.add(LSTM(25,return_sequences=True))
    model.add(Dropout(0.2))
    model.add(LSTM(10,return_sequences=True))
    model.add(Dropout(0.2))
    model.add(TimeDistributed(Dense(1)))
    opt = keras.optimizers.RMSprop(lr=learning_rate) #better performance
    #opt = keras.optimizers.Adam(lr=0.0001) # not good performance
    model.compile(loss="mse", optimizer=opt, metrics=["mse"])
    model.summary()
    return model

#Many to many model. The length of the input is different from the output length    
def buildManyToManyModel2(shape, steps_after):
    model = Sequential()
    model.add(LSTM(input_dim=shape[2], input_length=shape[1],output_dim=25, return_sequences=False))
    model.add(RepeatVector(steps_after))
    model.add(LSTM(output_dim=10, return_sequences=True))
    model.add(TimeDistributed(Dense(1)))
    model.add(Activation('linear'))  
    #model.compile(loss="mse", optimizer="adam", metrics=["mean_squared_error"])
    opt = keras.optimizers.RMSprop(lr=0.0001) #better performance  
    model.compile(loss="mse", optimizer=opt, metrics=["mse"])
    
    return model
    
def simpleAvg(X, n_pre, n_post):
    Y = []
    x = X.tolist()
   
    for i in range(n_post):
        m = np.mean(x[i:i+n_pre])
        x.append(m)
        Y.append(m)
  
    return np.array(Y)

def predictBasedSimpleAvg(X, target_column_index, n_pre, n_post):
    y_predict = []
    x = X[:, :,target_column_index]
    
    for i in range(X.shape[0]):
        y = simpleAvg(x[i], n_pre, n_post)
        y_predict.append(y)
    
    return np.array(y_predict)
        
def predict(model, input):
    y_pre = []
    for i in range(input.shape[0]):
        y_pre_i = model.predict(input[i])
        y_pre.append(y_pre_i)
    return np.array(y_pre)

def evaluateWithMetrics(true, lstm_predict, reg_predict, simple_avg_predict):
    sum_lstm_evs = 0
    sum_reg_evs = 0
    sum_simple_avg_evs = 0
    sum_lstm_mse = 0
    sum_reg_mse = 0
    sum_simple_avg_mse = 0
    sum_lstm_mae = 0
    sum_reg_mae = 0
    sum_simple_avg_mae = 0
    sum_lstm_r2_score = 0
    sum_reg_r2_score = 0
    sum_simple_avg_r2_score = 0
    
    lstm_evs = []
    reg_evs = []
    avg_evs = []
    
    lstm_mse = []
    reg_mse = []
    avg_mse = []
    
    lstm_mae = []
    reg_mae = []
    avg_mae = []
    
    lstm_r2_score = []
    reg_r2_score = []
    avg_r2_score = []
       
    for i in range(true.shape[0]):
        r = explained_variance_score(true[i], lstm_predict[i])
        lstm_evs.append(r)
        sum_lstm_evs = sum_lstm_evs + r
        
        r = explained_variance_score(true[i], reg_predict[i])
        reg_evs.append(r)
        sum_reg_evs = sum_reg_evs + r
        
        r = explained_variance_score(true[i], simple_avg_predict[i])
        avg_evs.append(r)
        sum_simple_avg_evs = sum_simple_avg_evs + r
        
        r = mean_squared_error(true[i], lstm_predict[i])
        lstm_mse.append(r)
        sum_lstm_mse = sum_lstm_mse + r
        
        r = mean_squared_error(true[i], reg_predict[i])
        reg_mse.append(r)
        sum_reg_mse = sum_reg_mse + r
        
        r = mean_squared_error(true[i], simple_avg_predict[i])
        avg_mse.append(r)       
        sum_simple_avg_mse = sum_simple_avg_mse + r
        
        r = mean_absolute_error(true[i], lstm_predict[i])
        lstm_mae.append(r)
        sum_lstm_mae = sum_lstm_mae + r
        
        r = mean_absolute_error(true[i], reg_predict[i])
        reg_mae.append(r)
        sum_reg_mae = sum_reg_mae + r
        
        r = mean_absolute_error(true[i], simple_avg_predict[i])
        avg_mae.append(r)
        sum_simple_avg_mae = sum_simple_avg_mae + r
        
        r = r2_score(true[i], lstm_predict[i])
        lstm_r2_score.append(r)
        sum_lstm_r2_score = sum_lstm_r2_score + r
        
        r = r2_score(true[i], reg_predict[i])
        reg_r2_score.append(r)
        sum_reg_r2_score = sum_reg_r2_score + r
        
        r = r2_score(true[i], simple_avg_predict[i])
        avg_r2_score.append(r)
        sum_simple_avg_r2_score = sum_simple_avg_r2_score + r
        print(f'mae:lstm:{mean_absolute_error(true[i], lstm_predict[i])}, avg:{mean_absolute_error(true[i], simple_avg_predict[i])}')
     
    plotMetrics2(lstm_evs, reg_evs, avg_evs, lstm_mse, reg_mse, avg_mse, lstm_mae, reg_mae, avg_mae, lstm_r2_score, reg_r2_score, avg_r2_score)
    
    avg_lstm_evs = sum_lstm_evs / true.shape[0]
    avg_reg_evs = sum_reg_evs / true.shape[0]
    avg_simple_avg_evs = sum_simple_avg_evs / true.shape[0]
    avg_lstm_mse = sum_lstm_mse / true.shape[0]
    avg_reg_mse = sum_reg_mse / true.shape[0]
    avg_simple_avg_mse = sum_simple_avg_mse / true.shape[0]
    avg_lstm_mae = sum_lstm_mae / true.shape[0]
    avg_reg_mae = sum_reg_mae / true.shape[0]
    avg_simple_avg_mae = sum_simple_avg_mae / true.shape[0]
    avg_lstm_r2_score = sum_lstm_r2_score / true.shape[0]
    avg_reg_r2_score = sum_reg_r2_score / true.shape[0]
    avg_simple_avg_r2_score = sum_simple_avg_r2_score / true.shape[0]
    
    print(f'explained variance score: lstm:{avg_lstm_evs}, regression:{avg_reg_evs}, simple avg:{avg_simple_avg_evs}')
    print(f'mean absolute error: lstm:{avg_lstm_mae}, regression:{avg_reg_mae}, simple avg: {avg_simple_avg_mae}')
    print(f'mean squared error: lstm:{avg_lstm_mse}, regression:{avg_reg_mse}, simple avg: {avg_simple_avg_mse}')
    print(f'r2 score: lstm:{avg_lstm_r2_score}, regression:{avg_reg_r2_score}, simple avg: {avg_simple_avg_r2_score}')
        

def mkdir(mypath):
    from errno import EEXIST
    from os import makedirs,path
    '''Creates a directory. equivalent to using mkdir -p on the command line'''
    try:
        makedirs(mypath)
    except OSError as exc:
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise
        
def plotManyToManyImage(X, Y, dates, predict, reg_predict, simple_avg_predict, target_index, num, folder_name):
    n_pre = X.shape[1]
    n_post = Y.shape[1]
    
    #Create new directory 
    mkdir(folder_name)
    
    nan_array = np.empty((n_pre - 1))
    nan_array.fill(np.nan)
    nan_array2 = np.empty(n_post)
    nan_array2.fill(np.nan)
    ind = np.arange(n_pre + n_post)
    fig, ax = plt.subplots()
       
    for i in range(num):
        forecasts = np.concatenate((nan_array, X[i, -1:, target_index], predict[i, :, 0]))
        simple_reg = np.concatenate((nan_array, X[i, -1:, target_index], reg_predict[i, :, 0]))
        ground_truth = np.concatenate((nan_array, X[i, -1:, target_index], Y[i, :, 0]))
        #The shape of X[i, :, target_index] is (n_pre,) 
        simple_mean = np.concatenate((nan_array, X[i, -1:, target_index], simple_avg_predict[i, :, 0]))
        network_input = np.concatenate((X[i, :, target_index], nan_array2))
        
        xs = [datetime.strptime(d, '%Y-%m-%d').date() for d in dates[i]]
        
        ax.plot(xs, network_input, 'b-x', label='Network input')
        ax.plot(xs, forecasts, 'r-x', label='LSTM forecast')
        ax.plot(xs, ground_truth, 'g-x', label = 'Ground truth')
        ax.plot(xs, simple_mean, 'm-x', label = 'Simple average forecast')
       # ax.plot(xs, simple_reg, 'y-x', label = 'Simple regression')
        #ax.xaxis.set_major_locator(MonthLocator())
        ax.grid(True)
        
        plt.xlabel('Dates')
        plt.ylabel('Exchange_rate')
        plt.title('Many to Many Forecast')
        plt.legend(loc='best')
        plt.gcf().autofmt_xdate()
        plt.savefig(folder_name + "/image" + '_' + str(i) + '.png')
        plt.cla()
        

def plotMetrics(lstm_evs, reg_evs, avg_evs, lstm_mse, reg_mse, avg_mse, lstm_mae, reg_mae, avg_mae, lstm_r2, reg_r2, avg_r2):
    
    plt.figure()
    plt.plot(lstm_evs, 'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_evs) for i in range(len(lstm_evs))])
    plt.plot(mean, 'r--', label = 'LSTM mean')
    plt.plot(reg_evs, 'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_evs) for i in range(len(reg_evs))])
    plt.plot(mean, 'b--', label = 'Simple Regression mean')
    plt.plot(avg_evs, 'g', label = 'Simple average')
    mean = np.array([np.mean(avg_evs) for i in range(len(avg_evs))])
    plt.plot(mean, 'g--', label = 'Simple average mean')
    
    plt.xlabel('Predict count')
    plt.ylabel('EVS')
    plt.title("Explained Variance Score")
    plt.legend(loc='best')
    plt.savefig("evs.png")
    
    plt.figure()
    plt.plot(lstm_mse,'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_mse) for i in range(len(lstm_mse))])
    plt.plot(mean,'r--', label = 'LSTM mean')
    plt.plot(reg_mse,'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_mse) for i in range(len(reg_mse))])
    plt.plot(mean,'b--', label = 'Simple Regression mean')
    plt.plot(avg_mse,'g', label = 'Simple average')
    mean = np.array([np.mean(avg_mse) for i in range(len(avg_mse))])
    plt.plot(mean,'g--', label = 'Simple average mean')
    plt.xlabel('Predict count')
    plt.ylabel('MSE')
    plt.title("Mean Squared Error")
    plt.legend(loc='best')
    plt.savefig("mse.png")
    
    plt.figure()
    plt.plot(lstm_mae,'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_mae) for i in range(len(lstm_mae))])
    plt.plot(mean,'r--', label = 'LSTM mean')
    plt.plot(reg_mae,'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_mae) for i in range(len(reg_mae))])
    plt.plot(mean,'b--', label = 'Simple Regression mean')
    plt.plot(avg_mae,'g', label = 'Simple average')
    mean = np.array([np.mean(avg_mae) for i in range(len(avg_mae))])
    plt.plot(mean,'g--', label = 'Simple average mean')
    plt.xlabel('Predict count')
    plt.ylabel('MAE')
    plt.title("Mean Absolute Error")
    plt.legend(loc='best')
    plt.savefig("mae.png")
    
    plt.figure()
    plt.plot(lstm_r2, 'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_r2) for i in range(len(lstm_r2))])
    plt.plot(mean, 'r--', label = 'LSTM mean')
    plt.plot(reg_r2, 'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_r2) for i in range(len(reg_r2))])
    plt.plot(mean, 'b--', label = 'Simple Regression mean')
    plt.plot(avg_r2, 'g',label = 'Simple average')
    mean = np.array([np.mean(avg_r2) for i in range(len(avg_r2))])
    plt.plot(mean, 'g--',label = 'Simple average mean')
    plt.xlabel('Predict count')
    plt.ylabel('R2 Score')
    plt.title("R Squared Score")
    plt.legend(loc='best')
    plt.savefig("r2.png")

def plotLoss(lstm, reg):
    plt.figure(figsize=(10, 4))
    plt.subplot(1,2,1)
    #plot loss change
    plt.plot(lstm.history['loss'], label = 'train loss')
    plt.plot(lstm.history['val_loss'], label = 'validation loss' )
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title("LSTM Model Loss")
    plt.legend(loc='best')
   
    plt.subplot(1,2,2)
    #plot loss change
    plt.plot(reg.history['loss'], label = 'train loss')
    plt.plot(reg.history['val_loss'], label = 'validation loss' )
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title("Simple Regression Model Loss")
    plt.legend(loc='best')
    plt.savefig("loss.png")
    

def plotLoss2(reg):
    plt.figure(figsize=(10, 4))
     
    #plot loss change
    plt.plot(reg.history['loss'], label = 'train loss')
    plt.plot(reg.history['val_loss'], label = 'validation loss' )
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title("Simple Regression Model Loss")
    plt.legend(loc='best')
    plt.show()
    
def plotMetrics2(lstm_evs, reg_evs, avg_evs, lstm_mse, reg_mse, avg_mse, lstm_mae, reg_mae, avg_mae, lstm_r2, reg_r2, avg_r2):
    
    plt.figure(figsize=(10, 10))
    plt.subplot(2,2,1)
    
    plt.plot(lstm_evs, 'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_evs) for i in range(len(lstm_evs))])
    plt.plot(mean, 'r--')
    plt.plot(reg_evs, 'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_evs) for i in range(len(reg_evs))])
    plt.plot(mean, 'b--')
    plt.plot(avg_evs, 'g', label = 'Simple average')
    mean = np.array([np.mean(avg_evs) for i in range(len(avg_evs))])
    plt.plot(mean, 'g--')
    plt.xlabel('Predict count')
    plt.ylabel('EVS')
    plt.title("Explained Variance Score")
    plt.legend(loc='best')
       
    plt.subplot(2,2,2)
    plt.plot(lstm_mse,'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_mse) for i in range(len(lstm_mse))])
    plt.plot(mean,'r--')
    plt.plot(reg_mse,'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_mse) for i in range(len(reg_mse))])
    plt.plot(mean,'b--')
    plt.plot(avg_mse,'g', label = 'Simple average')
    mean = np.array([np.mean(avg_mse) for i in range(len(avg_mse))])
    plt.plot(mean,'g--')
    plt.xlabel('Predict count')
    plt.ylabel('MSE')
    plt.title("Mean Squared Error")
    plt.legend(loc='best')
    
    plt.subplot(2,2,3)
    plt.plot(lstm_mae,'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_mae) for i in range(len(lstm_mae))])
    plt.plot(mean,'r--')
    plt.plot(reg_mae,'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_mae) for i in range(len(reg_mae))])
    plt.plot(mean,'b--')
    plt.plot(avg_mae,'g', label = 'Simple average')
    mean = np.array([np.mean(avg_mae) for i in range(len(avg_mae))])
    plt.plot(mean,'g--')
    plt.xlabel('Predict count')
    plt.ylabel('MAE')
    plt.title("Mean Absolute Error")
    plt.legend(loc='best')
    
    plt.subplot(2,2,4)
    plt.plot(lstm_r2, 'r', label = 'LSTM')
    mean = np.array([np.mean(lstm_r2) for i in range(len(lstm_r2))])
    plt.plot(mean, 'r--')
    plt.plot(reg_r2, 'b', label = 'Simple Regression')
    mean = np.array([np.mean(reg_r2) for i in range(len(reg_r2))])
    plt.plot(mean, 'b--')
    plt.plot(avg_r2, 'g',label = 'Simple average')
    mean = np.array([np.mean(avg_r2) for i in range(len(avg_r2))])
    plt.plot(mean, 'g--')
    plt.xlabel('Predict count')
    plt.ylabel('R2 Score')
    plt.title("R Squared Score")
    plt.legend(loc='best')
    plt.savefig("metrics.png")

#convert txt to csv  
#txt2csv("All_Data_Clean_1990-20170630_updated_11_27.txt", "exchange_rate.csv")

np.random.seed(7)

#read exchage_rate.csv
df = readData("exchange_rate.csv")
df_copy = df.copy()

#features need to be dropped
'''
drop_features = ['exchus', 'UK_3_Month_treasury_Securities', 'US_3_Month_Treasury_Yield', 'US_Federal_Funds_Effective', 'UK_London_Interbank_Offered_Rate',
                                                   'Interest_rate_required_reserves', 'inflation', 'totval', 'sprtrn', 'SP_500_risk_index', 
                                                   'efftax_Median', 'adv_sale_Median', 'staff_sale_Median', 'Employment_Cost_Index', 'Import_Price_Index',
                                                   'Export_Price_Index', 'Labor_Productivity_Index_business_sector', 'Labor_Productivity_Index_manufacturing_sector']

'''
drop_features = ['exchus', 'Labor_Productivity_Index_manufacturing_sector', 'sprtrn']

df = df.drop(drop_features, axis = 1)

#Normalization
df_norm = ZScoreNormalization(df)

#############################################################################################
# build Data, use last 12 months to predict next 1 month data
n_past = 48
n_future = 48
rate = 0.8
target_column_index, X_train, Y_train, dates = buildTrainData(df_norm, n_past, n_future)
# split training data and validation data
X_train, Y_train, X_val, Y_val, date_train, date_val = splitData(X_train, Y_train, dates, rate)
predicted_x, predicted_y, predicted_dates = buildOnePredictedData(df_norm, n_future, "2013-07")

#############################################################################################

reg_predicted_y = np.empty_like(Y_val)
reg_predicted_y.fill(np.nan)

'''

#############################################################################################
before_2008, after_2008 = chopDataSetByMonth(df_norm, "2008-02")
before_predict, predict_part = chopDataSetByMonth(df_norm, "2009-01")

#split time with 2008
target_column_index, X_train, Y_train, date_train = buildTrainData(before_2008, 12, 12)
target_column_index, X_val, Y_val, date_val = buildTrainData(predict_part, 12, 12)

predicted_x, predicted_y, predicted_dates = buildOnePredictedData(df_norm, 48, "2016-07")
#############################################################################################
'''
'''
#############################################################################################
before_2011, after_2011 = chopDataSetByMonth(df_norm, "2011-06")

#split time with 2008
target_column_index, X_train, Y_train, date_train = buildTrainData(before_2011, 48, 48)
predicted_x, predicted_y, predicted_dates = buildOnePredictedData(df_norm, 48, "2011-07")
#############################################################################################

'''

# from 2 dimmension to 3 dimension
Y_train = Y_train[:,:,np.newaxis]
Y_val = Y_val[:,:,np.newaxis]

'''
################################Grid Search########################################
epochs=[200, 300, 350, 400, 500, 1000]
batchs=[10, 20, 5]
learning_rates=[0.0001, 0.00001]

model = KerasRegressor(build_fn=buildManyToManyModel3, shape=X_train.shape)
param_grid = dict(epochs=epochs, batch_size=batchs, learning_rate=learning_rates)
grid = GridSearchCV(estimator=model, param_grid=param_grid)
grid_result = grid.fit(X_train, Y_train)

print(f"Best Run: {grid_result.best_score_} with {grid_result.best_params_}")
means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    print(f"{round(mean, 4)} ({round(stdev, 4)}) running with: {param}")
    
os._exit(0)
###################################################################################
'''

model = buildManyToManyModel3(X_train.shape)
#callback = EarlyStopping(monitor="loss", patience=10, verbose=1, mode="auto")
model.summary()

#fit model
#model.fit(X_train, Y_train, epochs=300, batch_size=128, validation_data=(X_val, Y_val), callbacks=[callback])
history = model.fit(X_train, Y_train, epochs=350, batch_size=10, validation_data=(X_val, Y_val), verbose = 2)

#predict with the model
y_predicted = model.predict(X_val)
future_y = model.predict(predicted_x)

#################################Build simple regression model###############################
#generate training data for simple regression model
simple_reg_train_data, simple_reg_val_data = chopDataSetByMonth(df_norm, date_val[0][n_future])
simple_reg_train_y = simple_reg_train_data['exukus'].values
simple_reg_train_x = simple_reg_train_data.drop(['Data', 'exukus'], axis = 1).values

#generate test data for simple regression model
sliding_win_size = n_future
simple_reg_val_x, simple_reg_val_y, simple_reg_dates = buildSimpleRegressionTestData(simple_reg_val_data, sliding_win_size)

#build a simple regression model
simple_reg_model = buildSimpleRegressionModel(simple_reg_train_x.shape[1])
reg_history = simple_reg_model.fit(simple_reg_train_x, simple_reg_train_y, epochs=600, batch_size=10, 
                                    validation_data=(simple_reg_val_x[0], simple_reg_val_y[0]), verbose = 2)

reg_predicted_y = predict(simple_reg_model, simple_reg_val_x)
#############################End Build simple regression model###############################


################################Predict using simple average#################################
simple_avg_predicted_y = predictBasedSimpleAvg(X_val, target_column_index, n_past, n_future)
simple_avg_predicted_y = simple_avg_predicted_y[:,:,np.newaxis]
######################################End simple average#####################################

'''
####################################Get layers Weights########################################## 
units = int(int(model.layers[0].trainable_weights[0].shape[1])/4)
print("No units: ", units)
W = model.layers[0].get_weights()[0]
U = model.layers[0].get_weights()[1]
b = model.layers[0].get_weights()[2]

W_i = W[:, :units]
W_f = W[:, units: units * 2]
W_c = W[:, units * 2: units * 3]
W_o = W[:, units * 3:]

print("W_i, W_f, W_c, W_o")
print(W_i, W_f, W_c, W_o)

print("U_i, U_f, U_c, U_o")

U_i = U[:, :units]
U_f = U[:, units: units * 2]
U_c = U[:, units * 2: units * 3]
U_o = U[:, units * 3:]
print(U_i, U_f, U_c, U_o)

print("b_i, b_f, b_c, b_o")
b_i = b[:units]
b_f = b[units: units * 2]
b_c = b[units * 2: units * 3]
b_o = b[units * 3:]
print(b_i, b_f, b_c, b_o)

###############################################################################################
'''
'''
###############################################################################################
folder_prefix = "use_2008_48_predict_maxmin"
###############################################################################################
'''

###############################################################################################
folder_prefix = "use_2008_48_predict_zscore_350"
###############################################################################################

save_path_past = folder_prefix + "_past"
save_path_future = folder_prefix + "_future"

plotManyToManyImage(X_val, Y_val, date_val, y_predicted, reg_predicted_y, simple_avg_predicted_y, target_column_index, Y_val.shape[0], save_path_past)
predicted_y = predicted_y[:,:, np.newaxis]
nan_array = np.empty_like(future_y)
nan_array.fill(np.nan)

plotManyToManyImage(predicted_x, predicted_y, predicted_dates, future_y, nan_array, nan_array, target_column_index, 1, save_path_future)
plotLoss(history, reg_history)
#plotLoss2(reg_history)

########################################Compare with metrics####################################
Y_val = np.squeeze(Y_val, axis=(2,))
y_predicted = np.squeeze(y_predicted, axis=(2,))
reg_predicted_y = np.squeeze(reg_predicted_y, axis=(2,))
simple_avg_predicted_y = np.squeeze(simple_avg_predicted_y, axis=(2,))

evaluateWithMetrics(Y_val, y_predicted, reg_predicted_y, simple_avg_predicted_y)
################################################################################################


