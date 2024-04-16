import os
import argparse
import pandas as pd
import joblib
import warnings

from sklearn.model_selection import *
from sklearn.metrics import *

warnings.filterwarnings("ignore", category=UserWarning, module='xgboost.core')


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
    if features and target_name:
        sel_features = pd.read_csv(features).iloc[:, 0].tolist()
        data = data[sel_features + [target_name]]

    return data


def split_data(data, target_name):
    '''
    Splits the dataset into training and testing sets.
    
    Args:
      data (DataFrame): The input dataset.
      target_name (str): Name of the target varivectorable column.
    
    Returns:
      tuple: Contains the training features (DataFrame), testing features (DataFrame), 
             training target (Series), and testing target (Series) datasets.
    '''
    features = data.drop(target_name, axis=1)
    target = data[target_name]
    features_train, features_test, target_train, target_test = train_test_split(features, target, test_size=0.1, random_state=42)

    return features_train, features_test, target_train, target_test



def train_model(model_untrained, features_train, target_train):
    '''
    Loads an untrained model from a file, trains it with the provided data, and returns the trained model.
    
    Args:
        model_untrained (str): Path to the saved untrained model file (.joblib).
        features_train (DataFrame): Training features.
        target_train (Series): Training target values.
        
    Returns:
        A trained model.
    '''
    
    model = joblib.load(model_untrained)
    model.fit(features_train, target_train)
        
    return model


def test_model(model, features_test, target_test):
    '''
    Tests a trained model on a test dataset to evaluate its performance.
    Calculates and returns the area under the Precision-Recall Curve (auPR) and
    the area under the Receiver Operating Characteristic Curve (auROC).

    Args:
        model: Trained model
        features_test (DataFrame): Test dataset features to make predictions on.
        target_test (Series): Actual target values for the test dataset.

    Returns:
        tuple: Contains the auPR and auROC scores (float).
    '''

    target_test_pred = model.predict(features_test)
    target_test_pred_probs = model.predict_proba(features_test)[:, 1]
    
    # calculate auPR
    precision, recall, _ = precision_recall_curve(target_test, target_test_pred_probs)
    auPR = auc(recall, precision)
    # calculate auROC
    auROC = roc_auc_score(target_test, target_test_pred)
    
    return auPR, auROC


if __name__ == "__main__":

    # define and parse command-line options:
    parser = argparse.ArgumentParser(description='Trains and evaluates a machine learning model using a specified dataset. \
                                                  The script facilitates the selection of specific features for training, loads an untrained model, trains it, \
                                                  and evaluates its performance on a test set. Additionally, it supports comparing against a baseline predictor if provided. \
                                                  Performance is evaluated using the area under the Precision-Recall Curve (auPR) and the area under the Receiver Operating Characteristic Curve (auROC).')

    parser.add_argument('--input', type=check_file_path, required=True, help='Path to the .csv file containing the dataset including features and target vector.')
    parser.add_argument('--features', type=check_file_path, required=True, help='Path to the .csv file listing the features to select, with features are listed in the first column.')
    parser.add_argument('--target', type=str, required=True, help='Name of the target vector column of input dataset.')
    parser.add_argument('--model', type=check_file_path, required=True, help='Path to the saved untrained model file (.joblib).')
    parser.add_argument('--output', default='trained_model.joblib', help='Path to output file. Defaults to "trained_model.joblib" in the current working directory.')
    parser.add_argument('--base', type=check_file_path, help='Optional: Path to .csv of baseline predictor dataset (also including target vector besides baseline feature).')

    args = parser.parse_args()

    # main
    data = load_select_features(args.input, args.features, args.target)
    features_train, features_test, target_train, target_test = split_data(data, args.target)

    model = train_model(args.model, features_train, target_train)
    joblib.dump(model, args.output)

    model_auPR, model_auROC = test_model(model, features_test, target_test)
    print('model - auPR:', model_auPR)
    print('model - auROC:', model_auROC)

    if args.base:
        base = load_select_features(args.base)
        base_train, base_test, target_train, target_test = split_data(base, args.target)

        base_model = train_model(args.model, base_train, target_train)

        base_auPR, base_auROC = test_model(base_model, base_test, target_test)
        print('baseline - auPR:', base_auPR)
        print('baseline - auROC:', base_auROC)






