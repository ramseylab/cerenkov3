import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils.multiclass import unique_labels
from sklearn.utils.extmath import softmax
from sklearn.utils.metaestimators import if_delegate_has_method


class Cerenkov3Classifier(BaseEstimator, ClassifierMixin):
    """ An example classifier which implements a 1-NN algorithm.
    For more information regarding how to build your own classifier, read more
    in the :ref:`User Guide <user_guide>`.
    Parameters
    ----------
    demo_param : str, default='demo'
        A parameter used for demonstation of how to pass and store paramters.
    Attributes
    ----------
    X_ : ndarray, shape (n_samples, n_features)
        The input passed during :meth:`fit`.
    y_ : ndarray, shape (n_samples,)
        The labels passed during :meth:`fit`.
    classes_ : ndarray, shape (n_classes,)
        The classes seen at :meth:`fit`.
    """

    def __init__(self, sgn_manager, n2v_manager, ext_manager, classifier, 
                 w_4DGp=1.0, w_4DGt=1.0, w_GTEx=1.0, w_TFBS=1.0, w_NG=1.0, w_coexp=1.0, w_hn=1.0, w_bg=1.0, 
                 d=4, l=10, r=10, k=5, p=1, q=1, w=True):
        self.sgn_manager = sgn_manager
        self.n2v_manager = n2v_manager
        self.ext_manager = ext_manager

        self.classifier = classifier

        self.w_4DGp = w_4DGp
        self.w_4DGt = w_4DGt
        self.w_GTEx = w_GTEx
        self.w_TFBS = w_TFBS
        self.w_NG = w_NG
        self.w_coexp = w_coexp
        self.w_hn = w_hn
        self.w_bg = w_bg

        self.d = d
        self.l = l
        self.r = r
        self.k = k
        self.p = p
        self.q = q
        self.w = w

    """
    DO NOT pass `ext_manager.load_feat()` as a `fit_param`, 
    because `fit_param` will be parititioned by CV
    """
    def fit(self, X_INT_ID, y):
        """A reference implementation of a fitting function for a classifier.
        Parameters
        ----------
        X_INT_ID : array-like, shape (n_samples, 1)
            The int IDs of training input samples.
        y : array-like, shape (n_samples,)
            The target values. An array of int.
        Returns
        -------
        self : object
            Returns self.
        """
        # Check that X and y have correct shape
        # X, y = check_X_y(X, y)
        # Store the classes seen during fit
        self.classes_ = unique_labels(y)

        # The order is [w_4DGp, w_4DGt, w_GTEx, w_TFBS, w_NG, w_coexp, w_hn, w_bg]
        network_path = self.sgn_manager.generate_sgn(w_4DGp=self.w_4DGp, 
                                                     w_4DGt=self.w_4DGt, 
                                                     w_GTEx=self.w_GTEx,
                                                     w_TFBS=self.w_TFBS, 
                                                     w_NG=self.w_NG,
                                                     w_coexp=self.w_coexp, 
                                                     w_hn=self.w_hn, 
                                                     w_bg=self.w_bg)
        n2v_feat_path = self.n2v_manager.generate_feat(network_path, 
                                                       d=self.d, 
                                                       l=self.l, 
                                                       r=self.r, 
                                                       k=self.k, 
                                                       p=self.p, 
                                                       q=self.q, 
                                                       w=self.w)
        # n2v_feat_df has columns `label_id_colname`, `int_id_colname` and `label_colname` (e.g. "name", "INT_ID" and "label")
        n2v_feat_df = pd.read_csv(n2v_feat_path, sep="\t")
        
        n2v_id_col = self.n2v_manager.label_id_colname
        n2v_label_col = self.n2v_manager.label_colname
        n2v_int_id_col = self.n2v_manager.int_id_colname

        if self.ext_manager:
            # ext_feat_df has columns `feat_id_colname`, `int_id_colname` and `feat_label_colname` (e.g. "name", "INT_ID" and "label")
            ext_feat_df = self.ext_manager.load_feat()

            ext_id_col = self.ext_manager.feat_id_colname
            ext_label_col = self.ext_manager.feat_label_colname
            ext_int_id_col = self.ext_manager.int_id_colname

            # trim ext_feat_df
            if ext_id_col in ext_feat_df:
                ext_feat_df.drop(columns=ext_id_col, inplace=True)
            if ext_label_col in ext_feat_df:
                ext_feat_df.drop(columns=ext_label_col, inplace=True)

            # merge ext_feat_df and n2v_feat_df on "INT_ID"
            if ext_int_id_col == n2v_int_id_col:
                self.feat_df_ = n2v_feat_df.merge(ext_feat_df, how="left", on=n2v_int_id_col)
            else:
                self.feat_df_ = n2v_feat_df.merge(ext_feat_df, how="left", left_on=n2v_int_id_col, right_on=ext_int_id_col)
                self.feat_df_.drop(columns=ext_int_id_col, inplace=True)
        else:
            self.feat_df_ = n2v_feat_df
            
        self.feat_df_.set_index(n2v_int_id_col, inplace=True)

        # dynamically determine the features after training starts
        # `self.X_` is the true feature matrix
        self.X_ = self.feat_df_.loc[X_INT_ID, :]

        if n2v_id_col in self.X_:
            self.X_.drop(columns=n2v_id_col, inplace=True)
        if n2v_label_col in self.X_:
            self.y_ = self.X_.loc[:, n2v_label_col]
            if not np.array_equal(self.y_, y):
                raise ValueError("input y differ from labels")
            self.X_.drop(columns=n2v_label_col, inplace=True)
        else:
            self.y_ = y
        
        # Return the classifier
        self.classifier.fit(self.X_, self.y_)

        return self

    def predict(self, X_INT_ID):
        """ A reference implementation of a prediction for a classifier.
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.
        Returns
        -------
        y : ndarray, shape (n_samples,)
            The label for each sample is the label of the closest sample
            seen during fit.
        """
        # Check is fit had been called
        check_is_fitted(self, ['X_', 'y_'])

        # Input validation
        # X = check_array(X)

        n2v_id_col = self.n2v_manager.label_id_colname
        n2v_label_col = self.n2v_manager.label_colname

        X_ = self.feat_df_.loc[X_INT_ID, :]
        if n2v_id_col in X_:
            X_.drop(columns=n2v_id_col, inplace=True)
        if n2v_label_col in X_:
            X_.drop(columns=n2v_label_col, inplace=True)

        y_ = self.classifier.predict(X_)
        return y_

    def predict_proba(self, X_INT_ID):
        # Check is fit had been called
        check_is_fitted(self, ['X_', 'y_'])

        # Input validation
        # X = check_array(X)

        n2v_id_col = self.n2v_manager.label_id_colname
        n2v_label_col = self.n2v_manager.label_colname

        X_ = self.feat_df_.loc[X_INT_ID, :]
        if n2v_id_col in X_:
            X_.drop(columns=n2v_id_col, inplace=True)
        if n2v_label_col in X_:
            X_.drop(columns=n2v_label_col, inplace=True)

        if hasattr(self.classifier, "predict_proba"):
            y_ = self.classifier.predict_proba(X_)
            return y_
        else:
            decision = self.decision_function(X_)
            if decision.ndim == 1:
                # Workaround for multi_class="multinomial" and binary outcomes
                # which requires softmax prediction with only a 1D decision.
                decision_2d = np.c_[-decision, decision]
            else:
                decision_2d = decision
            return softmax(decision_2d, copy=False)
