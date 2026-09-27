# Adapted from MedCAI/AttriMIL (https://github.com/MedCAI/AttriMIL), Apache-2.0.
# Modified to read cohort-specific UNI feature H5 files and resolve LUAD/LUSC
# from the actual feature path; see ATTRIMIL_LICENSE.
import os
import torch
import numpy as np
import pandas as pd
import math
import re
import pdb
import pickle
from scipy import stats

from torch.utils.data import Dataset
import h5py

def save_splits(split_datasets, column_keys, filename, boolean_style=False):
    splits = [split_datasets[i].slide_data['slide_id'] for i in range(len(split_datasets))]
    if not boolean_style:
        df = pd.concat(splits, ignore_index=True, axis=1)
        df.columns = column_keys
    else:
        df = pd.concat(splits, ignore_index = True, axis=0)
        index = df.values.tolist()
        one_hot = np.eye(len(split_datasets)).astype(bool)
        bool_array = np.repeat(one_hot, [len(dset) for dset in split_datasets], axis=0)
        df = pd.DataFrame(bool_array, index=index, columns = ['train', 'val', 'test'])

    df.to_csv(filename)
    print()

class Generic_WSI_Classification_Dataset(Dataset):
    def __init__(self,
        csv_path = 'dataset_csv/ccrcc_clean.csv',
        shuffle = False,
        seed = 7,
        print_info = True,
        label_dict = {},
        filter_dict = {},
        ignore=[],
        patient_strat=False,
        label_col = None,
        patient_voting = 'max',
        ):
        self.label_dict = label_dict
        self.num_classes = len(set(self.label_dict.values()))
        self.seed = seed
        self.print_info = print_info
        self.patient_strat = patient_strat
        self.train_ids, self.val_ids, self.test_ids  = (None, None, None)
        self.data_dir = None
        if not label_col:
            label_col = 'label'
        self.label_col = label_col

        slide_data = pd.read_csv(csv_path)
        slide_data = self.filter_df(slide_data, filter_dict)
        slide_data = self.df_prep(slide_data, self.label_dict, ignore, self.label_col)

        if shuffle:
            np.random.seed(seed)
            np.random.shuffle(slide_data)

        self.slide_data = slide_data
        self.patient_data_prep(patient_voting)
        self.cls_ids_prep()

        if print_info:
            self.summarize()

    def cls_ids_prep(self):
        self.patient_cls_ids = [[] for i in range(self.num_classes)]
        for i in range(self.num_classes):
            self.patient_cls_ids[i] = np.where(self.patient_data['label'] == i)[0]

        self.slide_cls_ids = [[] for i in range(self.num_classes)]
        for i in range(self.num_classes):
            self.slide_cls_ids[i] = np.where(self.slide_data['label'] == i)[0]

    def patient_data_prep(self, patient_voting='max'):
        patients = np.unique(np.array(self.slide_data['case_id']))
        patient_labels = []
        for p in patients:
            locations = self.slide_data[self.slide_data['case_id'] == p].index.tolist()
            assert len(locations) > 0
            label = self.slide_data['label'][locations].values
            if patient_voting == 'max':
                label = label.max()
            elif patient_voting == 'maj':
                label = stats.mode(label)[0]
            else:
                raise NotImplementedError
            patient_labels.append(label)
        self.patient_data = {'case_id':patients, 'label':np.array(patient_labels)}

    @staticmethod
    def df_prep(data, label_dict, ignore, label_col):
        if label_col != 'label':
            data['label'] = data[label_col].copy()
        mask = data['label'].isin(ignore)
        data = data[~mask]
        data.reset_index(drop=True, inplace=True)
        for i in data.index:
            key = data.loc[i, 'label']
            data.at[i, 'label'] = label_dict[key]
        return data

    def filter_df(self, df, filter_dict={}):
        if len(filter_dict) > 0:
            filter_mask = np.full(len(df), True, bool)
            for key, val in filter_dict.items():
                mask = df[key].isin(val)
                filter_mask = np.logical_and(filter_mask, mask)
            df = df[filter_mask]
        return df

    def __len__(self):
        if self.patient_strat:
            return len(self.patient_data['case_id'])
        else:
            return len(self.slide_data)

    def summarize(self):
        print("label column: {}".format(self.label_col))
        print("label dictionary: {}".format(self.label_dict))
        print("number of classes: {}".format(self.num_classes))
        for i in range(self.num_classes):
            print('Patient-LVL; Number of samples registered in class %d: %d' % (i, self.patient_cls_ids[i].shape[0]))
            print('Slide-LVL; Number of samples registered in class %d: %d' % (i, self.slide_cls_ids[i].shape[0]))

    def get_split_from_df(self, all_splits, split_key='train'):
        split = all_splits[split_key]
        split = split.dropna().reset_index(drop=True)
        if len(split) > 0:
            mask = self.slide_data['slide_id'].isin(split.tolist())
            df_slice = self.slide_data[mask].reset_index(drop=True)
            split = Generic_Split(df_slice, data_dir=self.data_dir, num_classes=self.num_classes)
        else:
            split = None
        return split

    def return_splits(self, from_id=True, csv_path=None):
        if from_id:
            train_split = Generic_Split(self.slide_data.loc[self.train_ids].reset_index(drop=True), self.data_dir, self.num_classes) if len(self.train_ids) > 0 else None
            val_split = Generic_Split(self.slide_data.loc[self.val_ids].reset_index(drop=True), self.data_dir, self.num_classes) if len(self.val_ids) > 0 else None
            test_split = Generic_Split(self.slide_data.loc[self.test_ids].reset_index(drop=True), self.data_dir, self.num_classes) if len(self.test_ids) > 0 else None
        else:
            assert csv_path
            all_splits = pd.read_csv(csv_path, dtype=self.slide_data['slide_id'].dtype)
            train_split = self.get_split_from_df(all_splits, 'train')
            val_split = self.get_split_from_df(all_splits, 'val')
            test_split = self.get_split_from_df(all_splits, 'test')
        return train_split, val_split, test_split

    def get_nearest(self, coords, k=8):
        from sklearn.neighbors import NearestNeighbors
        if torch.is_tensor(coords):
            coords = coords.cpu().numpy()
        nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        return indices[:, 1:]

    def get_feature_path(self, idx):
        """Return the actual cohort and feature file for a slide."""
        slide_id = str(self.slide_data['slide_id'].iloc[idx]).strip()
        filename = slide_id if slide_id.endswith('.h5') else f'{slide_id}.h5'
        if not isinstance(self.data_dir, dict):
            path = os.path.join(self.data_dir, 'features_uni_v1', filename)
            if not os.path.isfile(path):
                raise FileNotFoundError(path)
            return None, path

        preferred = None
        for column in ('Patho', 'source', 'project'):
            if column in self.slide_data.columns:
                value = str(self.slide_data[column].iloc[idx]).upper()
                preferred = next((key for key in self.data_dir if key.upper() in value), None)
                if preferred is not None:
                    break
        matches = [(key, os.path.join(root, 'features_uni_v1', filename))
                   for key, root in self.data_dir.items()
                   if os.path.isfile(os.path.join(root, 'features_uni_v1', filename))]
        if preferred is not None:
            for key, path in matches:
                if key == preferred:
                    return key, path
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise FileNotFoundError(f'No feature H5 found for {slide_id}')
        raise ValueError(f'{slide_id} is present in multiple cohorts; provide Patho/source/project in the metadata CSV')

    # ========================== 核心修改部分 ==========================
    def __getitem__(self, idx):
        slide_id = str(self.slide_data['slide_id'][idx]).strip()
        label = self.slide_data['label'][idx]

        _, full_path = self.get_feature_path(idx)

        with h5py.File(full_path, 'r') as hdf5_file:
            features = hdf5_file['features'][:]
            coords = hdf5_file['coords'][:]

        nearest = self.get_nearest(coords)

        return torch.from_numpy(features).float(), \
               label, \
               torch.from_numpy(coords).float(), \
               torch.from_numpy(nearest).long()
    # =================================================================

class Generic_MIL_Dataset(Generic_WSI_Classification_Dataset):
    def __init__(self, data_dir, **kwargs):
        super(Generic_MIL_Dataset, self).__init__(**kwargs)
        self.data_dir = data_dir
        self.use_h5 = True

    def load_from_h5(self, toggle):
        self.use_h5 = toggle

    def __getitem__(self, idx):
        if not self.use_h5:
            slide_id = self.slide_data['slide_id'][idx]
            label = self.slide_data['label'][idx]

            # 同样应用动态 source 逻辑
            project = None
            for col in ['Patho', 'source', 'project']:
                if col in self.slide_data.columns:
                    project = self.slide_data[col].iloc[idx]
                    break

            if isinstance(self.data_dir, dict):
                source = project if project in self.data_dir else list(self.data_dir.keys())[0]
                data_dir = self.data_dir[source]
            else:
                data_dir = self.data_dir

            full_path = os.path.join(data_dir, 'pt_files', '{}.pt'.format(slide_id))
            features = torch.load(full_path)
            return features, label

        return super(Generic_MIL_Dataset, self).__getitem__(idx)

class Generic_Split(Generic_MIL_Dataset):
    def __init__(self, slide_data, data_dir=None, num_classes=2):
        self.use_h5 = True
        self.slide_data = slide_data
        self.data_dir = data_dir
        self.num_classes = num_classes
        self.slide_cls_ids = [[] for i in range(self.num_classes)]
        for i in range(self.num_classes):
            self.slide_cls_ids[i] = np.where(self.slide_data['label'] == i)[0]

    def __len__(self):
        return len(self.slide_data)
