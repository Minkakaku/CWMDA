import argparse
import os
import time as tm
import pandas as pd
import numpy as np
import warnings

from torch.nn.modules.activation import Sigmoid
warnings.filterwarnings("ignore")

import torch
import torch.nn as nn
import torch.nn.functional as ft
from torch.distributions import constraints
from torch.nn.modules.linear import Linear
from utils.custom_mlp import MLP, Exp, SigmoidM
from utils.scdata_cached import mkdir_p, setup_data_loader, SingleCellCached, label2class_encoder, transform_class2label

import pyro
import pyro.distributions as dist
from pyro.contrib.examples.util import print_and_log
from pyro.infer import SVI, JitTrace_ELBO, JitTraceEnum_ELBO, Trace_ELBO, TraceEnum_ELBO, config_enumerate
from pyro.optim import Adam, ExponentialLR

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, matthews_corrcoef

class scClassifier3(nn.Module):
    """
    This class encapsulates the parameters (neural networks) and models & guides needed to train a
    semi-supervised variational auto-encoder on single cell datasets
    :param output_size:  size of the tensor representing the class label
    :param input_size: size of the tensor representing the cell
    :param z_dim: size of the tensor representing the latent random variable z
    :param hidden_layers: a tuple (or list) of MLP layers to be used in the neural networks
                          representing the parameters of the distributions in our model
    :param use_cude: use GPUs for faster training
    :param aux_loss_multiplier: the multiplier to use with the auxiliary loss
    """

    def __init__(self,
                 output_size = 10,
                 input_size = 2000,
                 z_dim = 50,
                 hidden_layers = (500,),
                 class2label = None,
                 config_enum = None,
                 use_cuda = False,
                 aux_loss_multiplier = None,
    ):
        super().__init__()

        # initialize the class with all arguments provided to the constructor
        self.output_size = output_size
        self.input_size = input_size
        self.z_dim = z_dim
        self.hidden_layers = hidden_layers
        self.class2label = class2label
        self.label2class = class2label
        self.allow_broadcast = config_enum == 'parallel'
        self.use_cuda = use_cuda
        self.aux_loss_multiplier = aux_loss_multiplier

        # define and instantiate the neural networks representing
        # the parameters of various distributions in the model
        self.setup_networks()

    def setup_networks(self):
        z_dim = self.z_dim
        hidden_sizes = self.hidden_layers

        # define the neural networks used later in the model and the guide.
        self.encoder_zy = MLP(
            [self.input_size] + hidden_sizes + [[z_dim, z_dim]],
            activation =  nn.Softplus,
            output_activation = [None, Exp],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.encoder_zy_y = MLP(
            [self.z_dim] + hidden_sizes + [self.output_size],
            activation = nn.Softplus,
            output_activation = nn.Softmax,
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.encoder_zy_z = MLP(
            [self.z_dim + self.output_size] + hidden_sizes + [[z_dim, z_dim]],
            activation =  nn.Softplus,
            output_activation = [None, Exp],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.encoder_ls = MLP(
            [self.input_size + self.output_size] + hidden_sizes + [[1,1]],
            activation = nn.Softplus,
            output_activation = [nn.Softplus, nn.Softplus],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.decoder_zy = MLP(
            [self.z_dim + self.output_size] + hidden_sizes + [[z_dim, z_dim]],
            activation =  nn.Softplus,
            output_activation = [None, Exp],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.decoder_concentrate = MLP(
            [self.z_dim] + hidden_sizes + [self.input_size],
            activation = nn.Softplus,
            output_activation = nn.Sigmoid,
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.cutoff = nn.Threshold(1.0e-9, 1.0e-9)
        self.softmax = nn.Softmax()

        # using GPUs for faster training of the networks
        if self.use_cuda:
            self.cuda()
    
    def model(self, xs, ys = None):
        """
        The model corresponds to the following generative process:
        p(z_theta) = normal(0, I)
        p(z_librarysize) = gamma(10000,1)
        p(y|x) = categorical(I/10.)
        p(theta|y,z_theta) = dirichlet(concentrate(y,z_theta))
        p(l|z_librarysize) = poisson(z_librarysize)
        p(x|theta, l) = multinomial(theta, l)
        concentrate is given by a neural network `decoder`

        :param xs: a batch of vectors of gene counts from a cell
        :param ys: (optional) a batch of the class labels
        :return: None
        """
        # register this pytorch module and all of its sub-modules with pyro
        pyro.module('scc', self)

        batch_size = xs.size(0)
        options = dict(dtype = xs.dtype, device = xs.device)
        with pyro.plate('data'):
            # if the label y is supervised, sample from the constant prior, 
            # otherwise, observe the value
            alpha_prior = torch.ones(batch_size, self.output_size, **options) / (
                1.9 * self.output_size
            )
            zy_ys = pyro.sample('zy_y', dist.OneHotCategorical(alpha_prior), obs = ys)

            # sample z from the constant prior distribution
            prior_loc = torch.zeros(batch_size, self.z_dim, **options)
            prior_scale = torch.ones(batch_size, self.z_dim, **options)
            zy_zs = pyro.sample('zy_z', dist.Normal(prior_loc, prior_scale).to_event(1))

            zy_loc, zy_scale = self.decoder_zy([zy_zs, zy_ys])
            zys = pyro.sample('zy', dist.Normal(zy_loc, zy_scale).to_event(1))

            # sample library size effect l from beta distribution
            ls_loc = torch.ones(batch_size, **options)
            ls_scale = torch.ones(batch_size, **options)
            ls = pyro.sample('ls', dist.Weibull(ls_loc, ls_scale))
            ls = ls.unsqueeze(-1)
            
            # sample multinomial parameter theta from dirichlet distribution
            concentrate = self.decoder_concentrate(zys)
            concentrate = concentrate * ls
            concentrate = self.cutoff(concentrate)

            # finally, score the observation
            max_count = torch.ceil(xs.sum(1).max()).int().item()
            pyro.sample('x', dist.DirichletMultinomial(total_count = max_count, concentration = concentrate), obs = xs)

    def guide(self, xs, ys = None):
        """
        The guide corresponds to the following:
        q(y|x) = categorical(alpha(x))
        q(z_theta|x,y) = normal(loc_theta(x,y), scale_theta(x,y))
        q(z_librarysize|x) = Gamma(concentrate_librarysize(x), rate_librarysize(x))
        alpha is given by a neural network `encoder_zy_y`
        loc_theta, scale_theta is given by a neural network `encoder_zy_z`
        concentrate_librarysize, rate_librarysize is given by a neural network `encoder_librarysize`

        :param xs: a batch of vectors of gene counts from a cell
        :param ys: (optional) a batch of the class labels
        :return: None
        """
        # inform Pyro that the variables in the batch of xs, ys are conditionally independent
        batch_size = xs.size(0)
        with pyro.plate('data'):
            zy_loc, zy_scale = self.encoder_zy(xs)
            zys = pyro.sample('zy', dist.Normal(zy_loc, zy_scale).to_event(1))

            # if the class label (the digit) is not supervised, sample 
            # the digit with the variational distribution
            # q(y|x) = categorical(alpha(x))
            alpha = self.encoder_zy_y(zys)
            if ys is None:
                ys = pyro.sample('zy_y', dist.OneHotCategorical(alpha))

            # sample the latent random variable with the variational
            # distribution q(z_theta|x,y) = normal(loc_theta(x,y), scale_theta(x,y))
            z_loc, z_scale = self.encoder_zy_z([zys, ys])
            zs = pyro.sample('zy_z', dist.Normal(z_loc, z_scale).to_event(1))

            # sample the latent random variable with the variational
            # distribution
            ls_loc, ls_scale = self.encoder_ls([xs,ys])
            ls = pyro.sample('ls', dist.Weibull(ls_loc.squeeze(), ls_scale.squeeze()))

    def classifier(self, xs):
        """
        classify a cell (or a batch of cells)

        :param xs: a batch of vectors of gene counts from a cell
        :return: a batch of the corresponding class labels (as one-hots)
                 along with the class probabilities
        """
        # use the trained model q(y|x) = categorical(alpha(x))
        # compute all class probabilities for the cell(s)
        zys,_ = self.encoder_zy(xs)
        alpha = self.encoder_zy_y(zys)

        # get the index (digit) that corresponds to
        # the maximum predicted class probability
        res, ind = torch.topk(alpha, 1)

        # convert the digit(s) to one-hot tensor(s)
        ys = torch.zeros_like(alpha).scatter_(1, ind, 1.0)

        return ys

    def convert_to_label(self, ys):
        if ys.is_cuda:
            ys = transform_class2label(ys.cpu().detach().numpy(), self.class2label)
        else:
            ys = transform_class2label(ys.numpy(), self.class2label)
        return ys
    
    def classifier_with_probability(self, xs, use_latent = False):
        """
        classify a cell (or a batch of cells)

        :param xs: a batch of vectors of gene counts from a cell
        :return: a batch of the corresponding class labels (as one-hots)
                 along with the class probabilities
        """
        # use the trained model q(y|x) = categorical(alpha(x))
        # compute all class probabilities for the cell(s)
        zys,_ = self.encoder_zy(xs)
        alpha = self.encoder_zy_y(zys)

        # get the index (digit) that corresponds to
        # the maximum predicted class probability
        res, ind = torch.topk(alpha, 1)

        # convert the digit(s) to one-hot tensor(s)
        ys = torch.zeros_like(alpha).scatter_(1, ind, 1.0)       

        return ys, alpha

    def classifier_BOG(self, xs, C = 100):
        pass

    def latent_embedding(self, xs):
        """
        compute the latent embedding of a cell (or a batch of cells)

        :param xs: a batch of vectors of gene counts from a cell
        :return: a batch of the latent embeddings
        """
        zs,_ = self.encoder_zy(xs)

        return zs

    def denoised_expression(self, xs):
        """
        compute the denoised expression of a cell (or a batch of cells)

        :param xs: a batch of vectors of gene counts from a cell
        :param ys: (optional) a batch of the class labels
        :return: a batch of the latent embeddings
        """
        zs,_ = self.encoder_zy(xs)
        alpha = self.decoder_concentrate(zs)

        return alpha

    def model_classify(self, xs, ys = None):
        """
        this model is used to add auxiliary (supervised) loss as described in the
        Kingma et al., "Semi-Supervised Learning with Deep Generative Models".
        """
        # register all pytorch (sub)modules with pyro
        pyro.module('scc', self)

        # inform pyro that the variables in the batch of xs, ys are conditionally independent
        with pyro.plate('data'):
            # this here is the extra term to yield an auxiliary loss that we do gradient descent on
            if ys is not None:
                zs,_ = self.encoder_zy(xs)
                alpha = self.encoder_zy_y(zs)
                with pyro.poutine.scale(scale = 50 * self.aux_loss_multiplier):
                    ys_aux = pyro.sample('y_aux', dist.OneHotCategorical(alpha), obs = ys)

    def guide_classify(self, xs, ys = None):
        """
        dummy guide function to accompany model_classify in inference
        """
        pass

def run_inference_for_epoch(sup_data_loader, unsup_data_loader, losses, use_cuda=True):
    """
    runs the inference algorithm for an epoch
    returns the values of all losses separately on supervised and unsupervised parts
    """
    num_losses = len(losses)

    # compute number of batches for an epoch
    sup_batches = len(sup_data_loader)
    unsup_batches = len(unsup_data_loader) if unsup_data_loader is not None else 0

    # initialize variables to store loss values
    epoch_losses_sup = [0.0] * num_losses
    epoch_losses_unsup = [0.0] * num_losses

    # setup the iterators for training data loaders
    sup_iter = iter(sup_data_loader)
    unsup_iter = iter(unsup_data_loader) if unsup_data_loader is not None else None

    # supervised data
    for i in range(sup_batches):
        # extract the corresponding batch
        (xs, ys) = next(sup_iter)

        if use_cuda:
            xs = xs.cuda()
            #if not isinstance(ys, tuple):
            if isinstance(ys, torch.Tensor):
                ys = ys.cuda() 

        # run the inference for each loss with supervised data as arguments
        for loss_id in range(num_losses):
            new_loss = losses[loss_id].step(xs, ys)
            epoch_losses_sup[loss_id] += new_loss

    # unsupervised data
    if unsup_data_loader is not None:
        for i in range(unsup_batches):
            # extract the corresponding batch
            (xs, ys) = next(unsup_iter)

            if use_cuda:
                xs = xs.cuda()
                #if not isinstance(ys, tuple):
                if isinstance(ys, torch.Tensor):
                    ys = ys.cuda() 

            # run the inference for each loss with unsupervised data as arguments
            for loss_id in range(num_losses):
                new_loss = losses[loss_id].step(xs)
                epoch_losses_unsup[loss_id] += new_loss

    # return the values of all losses
    return epoch_losses_sup, epoch_losses_unsup


def get_accuracy(data_loader, classifier_fn, batch_size):
    """
    compute the accuracy over the supervised training set or the testing set
    """
    predictions, actuals = [], []

    # use the appropriate data loader
    for (xs, ys) in data_loader:
        # use classification function to compute all predictions for each batch
        predictions.append(classifier_fn(xs))
        actuals.append(ys)

    # compute the number of accurate predictions
    predictions = torch.cat(predictions, dim=0)
    actuals = torch.cat(actuals, dim=0)
    _, y = torch.topk(actuals, 1)
    _, yhat = torch.topk(predictions, 1)
    y = y.detach().cpu().numpy()
    yhat = yhat.detach().cpu().numpy()
    #accurate_preds = 0
    #for pred, act in zip(predictions, actuals):
    #    for i in range(pred.size(0)):
    #        v = torch.sum(pred[i] == act[i])
    #        accurate_preds += v.item() == len(pred[i])
    #
    # calculate the accuracy between 0 and 1
    #accuracy = (accurate_preds * 1.0) / (len(predictions) * batch_size)
    accuracy = accuracy_score(y, yhat)
    f1_macro = f1_score(y, yhat, average='macro')
    f1_weighted = f1_score(y, yhat, average='weighted')
    precision = precision_score(y, yhat, average='macro')
    recall = recall_score(y, yhat, average='macro')
    mcc = matthews_corrcoef(y, yhat)

    return accuracy, f1_macro, f1_weighted, precision, recall, mcc

def label2class_map(sup_label_file, unsup_label_file = None, use_unsup = False):
    sup_labels = pd.read_csv(sup_label_file, header=None).squeeze().values
    if unsup_label_file is not None:
        unsup_labels = pd.read_csv(unsup_label_file, header=None).squeeze().values
        all_labels = np.concatenate((sup_labels, unsup_labels))
    else:
        all_labels = sup_labels
        #if use_unsup:
        #    all_labels = np.concatenate((all_labels, ['Undetermined']))
    return label2class_encoder(all_labels)

def main(args):
    """
    run inference for scClassifier3

    :param args: arguments for scClassifier3
    :return: None
    """
    if args.seed is not None:
        pyro.set_rng_seed(args.seed)

    if args.float64:
        torch.set_default_dtype(torch.float64)

    # create label to class mapping function
    use_unsup  = args.unsup_data_file is not None
    label2class = label2class_map(args.sup_label_file, args.unsup_label_file, use_unsup)

    # prepare dataloaders
    data_loaders = {'sup':None, 'unsup':None, 'valid':None}
    sup_num, unsup_num = 0, 0
    if args.sup_data_file is not None:
        data_loaders['sup'], data_loaders['valid'] = setup_data_loader(
            SingleCellCached, args.sup_data_file, args.sup_label_file, label2class, 'sup', args.validation_fold, args.log_transform, args.cuda, args.float64, args.batch_size
        )
        sup_num = len(data_loaders['sup'])
    if args.unsup_data_file is not None:
        data_loaders['unsup'], _ = setup_data_loader(
            SingleCellCached, args.unsup_data_file, args.unsup_label_file, label2class, 'unsup', 0, args.log_transform, args.cuda, args.float64, args.batch_size
        )
        unsup_num = len(data_loaders['unsup'])

    output_size = data_loaders['sup'].dataset.dataset.num_classes
    input_size = data_loaders['sup'].dataset.dataset.data.shape[1]

    # batch_size: number of cells (and labels) to be considered in a batch
    scc = scClassifier3(
        output_size=output_size,
        input_size=input_size,
        z_dim=args.z_dim,
        hidden_layers=args.hidden_layers,
        class2label=label2class,
        use_cuda=args.cuda,
        config_enum=args.enum_discrete,
        aux_loss_multiplier=args.aux_loss_multiplier,
    )

    # setup the optimizer
    adam_params = {'lr': args.learning_rate, 'betas':(args.beta_1, 0.999), 'weight_decay': 0.005}
    #optimizer = Adam(adam_params)
    optimizer = torch.optim.Adam
    decayRate = args.decay_rate
    scheduler = ExponentialLR({'optimizer': optimizer, 'optim_args': adam_params, 'gamma': decayRate})

    pyro.clear_param_store()

    # set up the loss(es) for inference, wrapping the guide in config_enumerate builds the loss as a sum
    # by enumerating each class label form the sampled discrete categorical distribution in the model
    guide = config_enumerate(scc.guide, args.enum_discrete, expand = True)
    Elbo = JitTraceEnum_ELBO if args.jit else TraceEnum_ELBO
    elbo = Elbo(max_plate_nesting=1, strict_enumeration_warning=False)
    loss_basic = SVI(scc.model, guide, scheduler, loss = elbo)

    # build a list of all losses considered
    losses = [loss_basic]

    # aux_loss: whether to use the auxiliary loss from NIPS 14 papers (Kingma et al)
    if args.aux_loss:
        elbo = JitTrace_ELBO() if args.jit else Trace_ELBO()
        loss_aux = SVI(scc.model_classify, scc.guide_classify, scheduler, loss = elbo)
        losses.append(loss_aux)

    try:
        # setup the logger if a filename is provided
        logger = open(args.logfile, 'w') if args.logfile else None

        # initializing local variables to maintain the best validation accuracy
        # seen across epochs over the supervised training set
        # and the corresponding testing set and the state of the networks
        best_valid_acc = 0.0
        best_valid_f1m = 0.0
        best_valid_mcc = 0.0

        asso_valid_f1m = 0.0
        asso_valid_f1w = 0.0
        asso_valid_pre = 0.0
        asso_valid_rec = 0.0
        asso_valid_mcc = 0.0
        
        asso_unsup_acc = 0.0
        asso_unsup_f1m = 0.0
        asso_unsup_f1w = 0.0
        asso_unsup_pre = 0.0
        asso_unsup_rec = 0.0
        asso_unsup_mcc = 0.0
        
        tr_start=tm.time()
        # run inference for a certain number of epochs
        for i in range(0, args.num_epochs):
            ep_tr_start = tm.time()

            # get the losses for an epoch
            epoch_losses_sup, epoch_losses_unsup = run_inference_for_epoch(
                data_loaders['sup'], data_loaders['unsup'], losses, args.cuda
            )

            # compute average epoch losses i.e. losses per example
            avg_epoch_losses_sup = map(lambda v: v / sup_num, epoch_losses_sup)
            avg_epoch_losses_unsup = map(lambda v: v / unsup_num, epoch_losses_unsup) if unsup_num > 0 else [0] * len(epoch_losses_unsup)
            avg_epoch_losses_sup = map(lambda v: "{:.4f}".format(v), avg_epoch_losses_sup)
            avg_epoch_losses_unsup = map(lambda v: "{:.4f}".format(v), avg_epoch_losses_unsup)

            # store the loss
            str_loss_sup = " ".join(map(str, avg_epoch_losses_sup))
            str_loss_unsup = " ".join(map(str, avg_epoch_losses_unsup))

            str_print = "{} epoch: avg losses {}".format(
                i+1, "{} {}".format(str_loss_sup, str_loss_unsup)
            )

            validation_accuracy, validation_f1_macro, validation_f1_weighted, validation_precision, validation_recall, validation_mcc = get_accuracy(
                data_loaders["valid"], scc.classifier, args.batch_size
            )

            str_print += " validation accuracy {:.4f}".format(validation_accuracy)
            str_print += " F1 {:.4f}(macro) {:.4f}(weighted)".format(validation_f1_macro, validation_f1_weighted)
            str_print += " precision {:.4f} recall {:.4f}".format(validation_precision, validation_recall)
            str_print += " mcc {:.4f}".format(validation_mcc)

            if (args.unsup_label_file is not None) and (args.unsup_data_file is not None):
                unsup_accuracy, unsup_f1_macro, unsup_f1_weighted, unsup_precision, unsup_recall, unsup_mcc = get_accuracy(
                    data_loaders['unsup'], scc.classifier, args.batch_size
                )            

            ep_tr_time = tm.time() - ep_tr_start
            str_print += " elapsed {:.4f} seconds".format(ep_tr_time)

            # update the best validation accuracy and the state of the parent 
            # module (including the networks)
            #if best_valid_acc <= validation_accuracy:
            if best_valid_mcc <= validation_mcc:
                do_update = False
                if best_valid_acc < validation_accuracy:
                    do_update = True
                elif best_valid_f1m < validation_f1_macro:
                    do_update = True

                if do_update:
                    best_valid_acc = validation_accuracy
                    best_valid_f1m = validation_f1_macro
                    best_valid_mcc = validation_mcc
                    
                    asso_valid_f1m = validation_f1_macro
                    asso_valid_f1w = validation_f1_weighted
                    asso_valid_pre = validation_precision
                    asso_valid_rec = validation_recall
                    asso_valid_mcc = validation_mcc

                    if (args.unsup_label_file is not None) and (args.unsup_data_file is not None):
                        asso_unsup_acc = unsup_accuracy
                        asso_unsup_f1m = unsup_f1_macro
                        asso_unsup_f1w = unsup_f1_weighted
                        asso_unsup_pre = unsup_precision
                        asso_unsup_rec = unsup_recall
                        asso_unsup_mcc = unsup_mcc

                    if args.save_model is not None:
                        if args.best_accuracy:
                            torch.save(scc, args.save_model)

            if i%args.decay_epochs == 0:
                scheduler.step() 

            if (i+1) == args.num_epochs:
                if args.save_model is not None:
                    if not args.best_accuracy:
                        torch.save(scc, args.save_model)

            print_and_log(logger, str_print)

        tr_time=tm.time()-tr_start
        if args.runtime:
            print('running time: {} secs'.format(tr_time))

        print_and_log(
            logger,
            "best validation accuracy {:.4f}".format(
                best_valid_acc
            ),
        )
        if (args.unsup_label_file is not None) and (args.unsup_data_file is not None):
            print_and_log(
                logger,
                "unsup accuracy: {:.4f} \nF1: {:.4f}(macro) {:.4f}(weighted) \nprecision {:.4f} recall {:.4f} \nmcc {:.4f}".format(
                    asso_unsup_acc, asso_unsup_f1m, asso_unsup_f1w, asso_unsup_pre, asso_unsup_rec, asso_unsup_mcc
                ),
            )
        else:
            print_and_log(
                logger,
                "F1: {:.4f}(macro) {:.4f}(weighted) \nprecision {:.4f} recall {:.4f} \nmcc {:.4f}".format(
                    asso_valid_f1m, asso_valid_f1w, asso_valid_pre, asso_valid_rec, asso_valid_mcc
                ),
            )
    finally:
        # close the logger file object if we opened it earlier
        if args.logfile:
            logger.close()

EXAMPLE_RUN = (
    "example run: python scClassifier3.py --seed 0 --cuda -n 50 --aux-loss -alm 46 -enum parallel "
    "--sup-data-file <sup_data_file> --sup-label-file <sup_label_file> "
    "--unsup-data-file [unsup_data_file] --unsup-label-file [sup_label_file] "
    "-zd 50 -hl 500 -lr 0.0001 -b1 0.99 -bs 100 -log ./tmp.log"
)

if __name__ == "__main__":
    #assert pyro.__version__.startswith("1.7.0")

    parser = argparse.ArgumentParser(description="scClassifier3\n{}".format(EXAMPLE_RUN))

    parser.add_argument(
        "--cuda", action="store_true", help="use GPU(s) to speed up training"
    )
    parser.add_argument(
        "--jit", action="store_true", help="use PyTorch jit to speed up training"
    )
    parser.add_argument(
        "-n", "--num-epochs", default=40, type=int, help="number of epochs to run"
    )
    parser.add_argument(
        "--aux-loss",
        action="store_true",
        help="whether to use the auxiliary loss from NIPS 14 paper "
        "(Kingma et al). It is not used by default ",
    )
    parser.add_argument(
        "-alm",
        "--aux-loss-multiplier",
        default=46,
        type=float,
        help="the multiplier to use with the auxiliary loss",
    )
    parser.add_argument(
        "-enum",
        "--enum-discrete",
        default="parallel",
        help="parallel, sequential or none. uses parallel enumeration by default",
    )
    parser.add_argument(
        "--sup-data-file",
        default=None,
        type=str,
        help="the data file of the supervised data",
    )
    parser.add_argument(
        "--sup-label-file",
        default=None,
        type=str,
        help="the label file of the supervised data",
    )
    parser.add_argument(
        "--unsup-data-file",
        default=None,
        type=str,
        help="the data file of the unsupervised data",
    )
    parser.add_argument(
        "--unsup-label-file",
        default=None,
        type=str,
        help="the label file of the unsupervised data",
    )
    parser.add_argument(
        "-64",
        "--float64",
        action="store_true",
        help="use double float precision",
    )
    parser.add_argument(
        "-lt",
        "--log-transform",
        action="store_true",
        help="run log-transform on count data",
    )
    parser.add_argument(
        "-vf",
        "--validation-fold",
        default=5,
        type=float,
        help="one of the folds of the supervised data for validation",
    )
    parser.add_argument(
        "-zd",
        "--z-dim",
        default=50,
        type=int,
        help="size of the tensor representing the latent variable z "
        "variable (handwriting style for our MNIST dataset)",
    )
    parser.add_argument(
        "-hl",
        "--hidden-layers",
        nargs="+",
        default=[500],
        type=int,
        help="a tuple (or list) of MLP layers to be used in the neural networks "
        "representing the parameters of the distributions in our model",
    )
    parser.add_argument(
        "-lr",
        "--learning-rate",
        default=0.0001,
        type=float,
        help="learning rate for Adam optimizer",
    )
    parser.add_argument(
        "-dr",
        "--decay-rate",
        default=0.9,
        type=float,
        help="decay rate for Adam optimizer",
    )
    parser.add_argument(
        "-de",
        "--decay-epochs",
        default=20,
        type=int,
        help="decay learning rate every #epochs",
    )
    parser.add_argument(
        "-b1",
        "--beta-1",
        default=0.95,
        type=float,
        help="beta-1 parameter for Adam optimizer",
    )
    parser.add_argument(
        "-bs",
        "--batch-size",
        default=1000,
        type=int,
        help="number of images (and labels) to be considered in a batch",
    )
    parser.add_argument(
        "-rt",
        "--runtime",
        action="store_true",
        help="print running time",
    )
    parser.add_argument(
        "-log",
        "--logfile",
        default="./tmp.log",
        type=str,
        help="filename for logging the outputs",
    )
    parser.add_argument(
        "--seed",
        default=None,
        type=int,
        help="seed for controlling randomness in this example",
    )
    parser.add_argument(
        "--save-model",
        default=None,
        type=str,
        help="path to save model for prediction",
    )
    parser.add_argument(
        "-ba",
        "--best-accuracy",
        action="store_true",
        help="save the model with best classification accuracy",
    )
    args = parser.parse_args()

    assert (
        (args.sup_data_file is not None) and (os.path.exists(args.sup_data_file))
    ), "sup_data_file must be provided"
    assert (
        (args.sup_label_file is not None) and (os.path.exists(args.sup_label_file))
    ), "sup_data_file must be provided"
    assert (
        (args.validation_fold >= 0)
    ), "fold of the supervised data used for validation should be greater than 0"
    if args.validation_fold > 0:
        args.validation_fold = 1. / args.validation_fold
    
    main(args)
