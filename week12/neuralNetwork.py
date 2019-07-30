import platform
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras import layers
from sklearn import datasets
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from time import time
from tensorflow.keras.callbacks import TensorBoard
from sklearn.preprocessing import MinMaxScaler
print(tf.__version__)
print("machine: ", platform.machine())
print("OS version: ", platform.version())
print("OS description: ", platform.uname())
print("Processor: ", platform.processor())

data = pd.read_csv("/Users/farhanir/PycharmProjects/qtw/week12/nn/data/HIGGS.csv")

data.describe()
data.head()
list(data)
data.columns = ['class','lepton pT', 'lepton eta', 'lepton phi', 'missing energy magnitude', 'missing energy phi', 'jet 1 pt', 'jet 1 eta', 'jet 1 phi', 'jet 1 b-tag', 'jet 2 pt', 'jet 2 eta', 'jet 2 phi', 'jet 2 b-tag', 'jet 3 pt', 'jet 3 eta', 'jet 3 phi', 'jet 3 b-tag', 'jet 4 pt', 'jet 4 eta', 'jet 4 phi', 'jet 4 b-tag', 'm_jj', 'm_jjj', 'm_lv', 'm_jlv', 'm_bb', 'm_wbb', 'm_wwbb']

list(data)

x = data.drop('class', axis=1)

y = data['class'].as_matrix()

scaler = MinMaxScaler(feature_range=(0, 1))
scaled_train = scaler.fit_transform(x)

# Print out the adjustment that the scaler applied to the total_earnings column of data
print("Note: median values were scaled by multiplying by {:.10f} and adding {:.6f}".format(scaler.scale_[12], scaler.min_[12]))
multiplied_by = scaler.scale_[12]
added = scaler.min_[12]

scaled_train_df = pd.DataFrame(scaled_train, columns=x.columns.values)
for i in scaled_train_df:
    scaled_train_df[i].hist()
    plt.show()

model = tf.keras.Sequential()
model.add(layers.Dense(300, activation='tanh'))
model.add(layers.Dense(300, activation='tanh'))
model.add(layers.Dense(300, activation='tanh'))
model.add(layers.Dense(300, activation='tanh'))
model.add(layers.Dense(1))


def auroc(y_true, y_pred):
    return tf.py_func(roc_auc_score, (y_true, y_pred), tf.double)


tb = TensorBoard(log_dir=f"logs\\{time()}")
model.compile(optimizer='sgd', loss='mean_squared_error', metrics=['accuracy', auroc])
history = model.fit(scaled_train_df.values, y, epochs=50, batch_size=100, callbacks=[tb])

# Plot training & validation accuracy values
plt.plot(history.history['acc'])
plt.plot(history.history['auroc'])
plt.title('ROC AUC Score against in 10 Epochs')
plt.ylabel('ROC AUC Score')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='upper left')
plt.show()

#
# # Plot training & validation loss values
# plt.plot(history.history['loss'])
# plt.plot(history.history['auroc'])
# plt.title('Model loss')
# plt.ylabel('Loss')
# plt.xlabel('Epoch')
# plt.legend(['Train', 'Test'], loc='upper left')
# plt.show()



### References:
# 1. https://arxiv.org/pdf/1402.4735.pdf
# 2. https://archive.ics.uci.edu/ml/datasets/HIGGS
# 3. https://stackoverflow.com/questions/3103178/how-to-get-the-system-info-with-python
# 4. https://cv-tricks.com/tensorflow-tutorial/keras/
# 5. https://medium.com/ibm-data-science-experience/markdown-for-jupyter-notebooks-cheatsheet-386c05aeebed
