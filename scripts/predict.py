import os
import argparse
import pandas as pd
import joblib
from sklearn.metrics import *

# author: Susanne Voigt


# functions 
def check_file_path(file_path):
    '''
    Checks if the specified file path exists.

    Args:
      file_path (str): The path to the file.

    Returns:
      str: The validated file path if it exists.

    Raises:
      argparse.ArgumentTypeError: If the file path does not exist.
    '''
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"{file_path} does not exist.")

    return file_path
    

def load_select_features(input, features=None, target_name=None):
    '''
    Loads a dataset from a .csv file and selects specific features along with a target column.

    Args:
    - input (str): Path to the input .csv file containing the dataset.
    - features (str, optional): Path to the .csv file listing the features to select in a DataFrame. Assumes features are listed in the first column.
    - target_name (str, optional): Name of the target vector column.

    Returns:
    - DataFrame containing only the selected features and the target column.
    '''
    data = pd.read_csv(input, index_col=0)
    if features:
        sel_features = pd.read_csv(features).iloc[:, 0].tolist()
        if target_name:
            data = data[sel_features + [target_name]]
        else:
            data = data[sel_features]

    return data


def predict(model_trained, data, features=None):
    '''
    Generates predictions using the provided model and data. 

    Args:
        model (str):  Path to the saved trained model to use for making predictions.
        data (DataFrame): The dataset to predict on.
        features (list, optional): A list of column names to be used for making predictions. 
                                   If not provided, predictions are made using all columns in `data`.

    Returns:
        np.ndarray: An array of predictions made by the model.
    '''
    model = joblib.load(model_trained)
    
    if features:
        sel_features = pd.read_csv(features).iloc[:, 0].tolist()
        pred = model.predict(data[sel_features])
    else:
        pred = model.predict(data)

    return pred


def evaluate_pred(model_trained, data, target_name, features=None):
    '''
    Evaluates predictions, if target vector is available,.
    Calculates and returns the area under the Precision-Recall Curve (auPR) and
    the area under the Receiver Operating Characteristic Curve (auROC).

    Args:
        model (str):  Path to the saved trained model to use for making predictions.
        data (DataFrame): The dataset including target vector.
        target_name (str): Name of the target vector column.
        features (list, optional): A list of column names to be used for making predictions. 
                                   If not provided, predictions are made using all columns in `data`.

    Returns:
        tuple: Contains the auPR and auROC scores (float).
    '''

    model = joblib.load(model_trained)
    if features:
        sel_features = pd.read_csv(features).iloc[:, 0].tolist()
        pred_probs = model.predict_proba(data[sel_features])[:, 1]
        pred = model.predict(data[sel_features])
    else:
        pred_probs = model.predict_proba(data.drop(target_name, axis=1))[:, 1]
        pred = model.predict(data.drop(target_name, axis=1))
    
    # calculate auPR
    precision, recall, _ = precision_recall_curve(data[target_name], pred_probs)
    auPR = auc(recall, precision)
    # calculate auROC
    auROC = roc_auc_score(data[target_name], pred)
    
    return auPR, auROC


if __name__ == "__main__":

    # define and parse command-line options:
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--input', type=check_file_path, required=True, help='Path to the .csv file containing the dataset including features and target vector.')
    parser.add_argument('--features', type=check_file_path, required=True, help='Path to the .csv file listing the features to select, with features are listed in the first column.')
    parser.add_argument('--model', type=check_file_path, required=True, help='Path to the saved trained model file (.joblib).')
    parser.add_argument('--output', default='data_pred.csv', help='Path to output file. Defaults to "data_pred.csv" in the current working directory.')
    parser.add_argument('--target', type=str, help='Optional: Name of the target vector column of input dataset.')

    args = parser.parse_args()

    # main
    data = load_select_features(args.input, args.features, args.target)
    pred = predict(args.model, data, args.features)
    data['pred'] = pred
    data.to_csv(args.output)

    if args.target:
        auPR, auROC = evaluate_pred(args.model, data, args.target, args.features)
        print('auPR:', auPR)
        print('auROC:', auROC)
