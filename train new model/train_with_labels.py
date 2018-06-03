from __future__ import print_function
import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping,ModelCheckpoint
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import metrics

data_augmentation = False
# num_predictions = 20
batch_size = 1024
num_classes = 3
epochs = 200
data_augmentation = False
# num_predictions = 20
model_name = 'keras_cnn_trained_model_shallow.h5'
# The data, shuffled and split between train and test sets:

def load_data_TF2(indel_list): # cell type specific  ## random samples for reactome is not enough, need borrow some from keggp
    import random
    import numpy as np
    xxdata_list = []
    yydata = []
    count_set = [0]
    count_setx = 0
    for i in indel_list:#len(h_tf_sc)):
        xdata = np.load('data_label/Nxdata_tf' + str(i) + '.npy')
        ydata = np.load('data_label/ydata_tf' + str(i) + '.npy')
        for k in range(len(ydata)):
            xxdata_list.append(xdata[k,:,:,:])
            yydata.append(ydata[k])
        count_setx = count_setx + len(ydata)
        count_set.append(count_setx)
        print (i,len(ydata))
    yydata_array = np.array(yydata)
    yydata_x = yydata_array.astype('int')
    print (np.array(xxdata_list).shape)
    return((np.array(xxdata_list),yydata_x,count_set))

if len(sys.argv) < 3:
    print ('No enough input files')
    sys.exit()
length_TF =int(sys.argv[1])
whole_data_TF = [i for i in range(length_TF)]
for test_indel in range(1,4):
    test_TF = [i for i in range (int(np.ceil((test_indel-1)*0.333*length_TF)),int(np.ceil(test_indel*0.333*length_TF)))]
    train_TF = [i for i in whole_data_TF if i not in test_TF]
    (x_train, y_train,count_set_train) = load_data_TF2(train_TF)
    (x_test, y_test,count_set_test) = load_data_TF2(test_TF)
    print(x_train.shape, 'x_train samples')
    print(x_test.shape, 'x_test samples')
    save_dir = os.path.join(os.getcwd(),str(test_indel)+'_saved_models_T_32-32-64-64-128-128-512_e'+str(epochs)+'_threesets_acc_new_unique_xy_new_rand')
    y_train = keras.utils.to_categorical(y_train, num_classes)
    y_test = keras.utils.to_categorical(y_test, num_classes)
    print(y_train.shape, 'y_train samples')
    print(y_test.shape, 'y_test samples')
    ############
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    ############
    model = Sequential()
    model.add(Conv2D(32, (3, 3), padding='same',
                     input_shape=x_train.shape[1:]))
    model.add(Activation('relu'))
    model.add(Conv2D(32, (3, 3)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))

    model.add(Conv2D(64, (3, 3), padding='same'))
    model.add(Activation('relu'))
    model.add(Conv2D(64, (3, 3)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))

    model.add(Conv2D(128, (3, 3), padding='same'))
    model.add(Activation('relu'))
    model.add(Conv2D(128, (3, 3)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(512))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(num_classes))
    model.add(Activation('softmax'))
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(optimizer=sgd,
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])

    early_stopping = keras.callbacks.EarlyStopping(monitor='val_acc', patience=50, verbose=0, mode='auto')
    checkpoint1 = ModelCheckpoint(filepath=save_dir + '/weights.{epoch:02d}-{val_loss:.2f}.hdf5', monitor='val_loss',
                                  verbose=1, save_best_only=False, save_weights_only=False, mode='auto', period=1)
    checkpoint2 = ModelCheckpoint(filepath=save_dir + '/weights.hdf5', monitor='val_acc', verbose=1,
                                  save_best_only=True, mode='auto', period=1)
    callbacks_list = [checkpoint2, early_stopping]
    if not data_augmentation:
        print('Not using data augmentation.')
        history = model.fit(x_train, y_train,
                  batch_size=batch_size,
                  epochs=epochs,validation_split=0.2,
                  shuffle=True, callbacks=callbacks_list)

    # Save model and weights

    model_path = os.path.join(save_dir, model_name)
    model.save(model_path)
    print('Saved trained model at %s ' % model_path)
    # Score trained model.
    scores = model.evaluate(x_test, y_test, verbose=1)
    print('Test loss:', scores[0])
    print('Test accuracy:', scores[1])
    y_predict = model.predict(x_test)
    np.save(save_dir+'/end_y_test.npy',y_test)
    np.save(save_dir+'/end_y_predict.npy',y_predict)

    plt.figure(figsize=(10, 6))
    plt.subplot(1,2,1)
    plt.plot(history.history['acc'])
    plt.plot(history.history['val_acc'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.grid()
    plt.legend(['train', 'val'], loc='upper left')
    plt.subplot(1,2,2)
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.grid()
    plt.savefig(save_dir+'/end_result.png')
    plt.figure(figsize=(10, 6))
    for i in range(3):
        y_test_x = [j[i] for j in y_test]
        y_predict_x = [j[i] for j in y_predict]
        fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
        plt.subplot(1, 3, i + 1)
        plt.plot(fpr, tpr)
        plt.grid()
        plt.plot([0, 1], [0, 1])
        plt.xlabel('FP')
        plt.ylabel('TP')
        plt.ylim([0, 1])
        plt.xlim([0, 1])
        auc = np.trapz(tpr, fpr)
        print('AUC:', auc)
        plt.title('label' + str(i) + ', AUC:' + str(auc))
    plt.savefig(save_dir + '/end_3labels.png')
    plt.figure(figsize=(10, 6))
    y_predict_x = [j[1] + j[2] for j in y_predict]
    y_test_x = [1 - j[0] for j in y_test]
    fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
    # Print ROC curve
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1])
    # Print AUC
    auc = np.trapz(tpr, fpr)
    print('AUC:', auc)
    plt.ylim([0, 1])
    plt.xlim([0, 1])
    plt.grid()
    plt.title('label 1+2, AUC:' + str(auc))
    plt.xlabel('FP')
    plt.ylabel('TP')
    plt.savefig(save_dir + '/end_1+2.png')
    #######################################
    plt.figure(figsize=(10, 6))
    y_predict1 = []
    y_test1 = []
    x = 2
    for i in range(int(len(y_predict) / 3)):
        y_predict1.append(y_predict[3 * i][x] - y_predict[3 * i + 1][x])
        y_predict1.append(-y_predict[3 * i][x] + y_predict[3 * i + 1][x])
        y_test1.append(y_test[3 * i][x])
        y_test1.append(y_test[3 * i + 1][x])
    fpr, tpr, thresholds = metrics.roc_curve(y_test1, y_predict1, pos_label=1)
    # Print ROC curve
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1])
    # Print AUC
    auc = np.trapz(tpr, fpr)
    print('AUC:', auc)
    plt.ylim([0, 1])
    plt.xlim([0, 1])
    plt.grid()
    plt.title('label 1 vs 2,direction diff, AUC:' + str(auc))
    plt.xlabel('FP')
    plt.ylabel('TP')
    plt.savefig(save_dir + '/end_1vs2.png')
    ###########################################################3
    model.load_weights(save_dir + '/weights.hdf5')
    scores = model.evaluate(x_test, y_test, verbose=1)
    print('Test loss:', scores[0])
    print('Test accuracy:', scores[1])
    y_predict = model.predict(x_test)
    np.save(save_dir+'/min_y_test.npy',y_test)
    np.save(save_dir+'/min_y_predict.npy',y_predict)

    plt.figure(figsize=(10, 6))
    plt.subplot(1,2,1)
    plt.plot(history.history['acc'])
    plt.plot(history.history['val_acc'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.grid()
    plt.legend(['train', 'val'], loc='upper left')
    plt.subplot(1,2,2)
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.grid()
    plt.savefig(save_dir+'/min_result.png')
    plt.figure(figsize=(10, 6))
    for i in range(3):
        y_test_x = [j[i] for j in y_test]
        y_predict_x = [j[i] for j in y_predict]
        fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
        plt.subplot(1, 3, i + 1)
        plt.plot(fpr, tpr)
        plt.grid()
        plt.plot([0, 1], [0, 1])
        plt.xlabel('FP')
        plt.ylabel('TP')
        plt.ylim([0, 1])
        plt.xlim([0, 1])
        auc = np.trapz(tpr, fpr)
        print('AUC:', auc)
        plt.title('label' + str(i) + ', AUC:' + str(auc))
    plt.savefig(save_dir + '/min_3labels.png')
    plt.figure(figsize=(10, 6))
    y_predict_x = [j[1] + j[2] for j in y_predict]
    y_test_x = [1 - j[0] for j in y_test]
    fpr, tpr, thresholds = metrics.roc_curve(y_test_x, y_predict_x, pos_label=1)
    # Print ROC curve
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1])
    # Print AUC
    auc = np.trapz(tpr, fpr)
    print('AUC:', auc)
    plt.ylim([0, 1])
    plt.xlim([0, 1])
    plt.grid()
    plt.title('label 1+2, AUC:' + str(auc))
    plt.xlabel('FP')
    plt.ylabel('TP')
    plt.savefig(save_dir + '/min_1+2.png')
    #######################################
    plt.figure(figsize=(10, 6))
    y_predict1 = []
    y_test1 = []
    x = 2
    for i in range(int(len(y_predict) / 3)):
        y_predict1.append(y_predict[3 * i][x] - y_predict[3 * i + 1][x])
        y_predict1.append(-y_predict[3 * i][x] + y_predict[3 * i + 1][x])
        y_test1.append(y_test[3 * i][x])
        y_test1.append(y_test[3 * i + 1][x])
    fpr, tpr, thresholds = metrics.roc_curve(y_test1, y_predict1, pos_label=1)
    # Print ROC curve
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1])
    # Print AUC
    auc = np.trapz(tpr, fpr)
    print('AUC:', auc)
    plt.ylim([0, 1])
    plt.xlim([0, 1])
    plt.grid()
    plt.title('label 1 vs 2,direction diff, AUC:' + str(auc))
    plt.xlabel('FP')
    plt.ylabel('TP')
    plt.savefig(save_dir + '/min_1vs2.png')
