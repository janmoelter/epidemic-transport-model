# ******************************************************************************
# This code shall be licensed under the GNU General Public License v3.0.
#
# (c) 2021-2022 , Jan Mölter
# ******************************************************************************

__all__ = ['epidemic_transport_model']

import os
import sys

import warnings

from copy import deepcopy

import random as rd

import enum

import numpy as np
import networkx as nx
import pandas as pd

import heapq

import scipy.linalg
import scipy.integrate
import scipy.sparse

import math


class epidemic_transport_model(object):

    class event(object):

        class ACTION(enum.IntEnum):
            UNSPECIFIED = 0
            MOVE = 1
            COMMUNITY_INFECTION = 2
            TRANSPORT_INFECTION = 3
            RECOVERY = 4
            IMMUNITY_LOSS = 5

        def __init__(self, time, subject, sender, action):

            self.time = time
            self.subject = subject
            self.sender = sender
            self.action = action

        def __lt__(self, other):
            if not (self.time == other.time):
                return self.time < other.time
            else:
                return self.time <= other.time and self.action <= other.action

        def submit(self, heap_queue, filter=None):

            if filter is None:
                filter = lambda e: True

            if filter(self):
                heapq.heappush(heap_queue, self)

    def __init__(self, transport_networks, community_size, community_network_degree, initial_site, initial_prevalence, epidemic_parameters, fractional_exponent=1, transport_network_interpolation_functional=(lambda t : 0), skip_community_network_generation=False):
        
        if not (isinstance(transport_networks, (nx.classes.graph.Graph, nx.classes.digraph.DiGraph)) or (isinstance(transport_networks, (list, tuple)) and len(transport_networks) == 2 and all([isinstance(_, (nx.classes.graph.Graph, nx.classes.digraph.DiGraph)) for _ in transport_networks]))):
            raise TypeError('`transport_networks` is expected to be of type nx.classes.digraph.Graph or nx.classes.digraph.DiGraph or of type list with length 2 and elements of that type.')

        if not (isinstance(community_size, int)):
            raise TypeError('`community_size` is expected to be of type int.')

        if not (isinstance(community_network_degree, (int, float)) or (isinstance(community_network_degree, list) and all([isinstance(_, int) for _ in community_network_degree]))):
            raise TypeError('`community_network_degree` is expected to be of type int or float or of type list with elements of type int.')
        else:
            if isinstance(community_network_degree, (int, float)):
                community_network_degree = [community_network_degree] * community_size
                if isinstance(community_network_degree, float):
                    skip_community_network_generation = True

        if not (isinstance(initial_site, int)):
            raise TypeError('`initial_site` is expected to be of type int.')

        if not (isinstance(initial_prevalence, (float, int))):
            raise TypeError('`initial_prevalence` is expected to be of type float.')

        if not (isinstance(fractional_exponent, (float, int))):
            raise TypeError('`fractional_exponent` is expected to be of type float.')

        if not (isinstance(epidemic_parameters, dict) and all([isinstance(_, (float, int)) for _ in epidemic_parameters.values()])):
            raise TypeError('`epidemic_parameters` is expected to be of type dict with entries of type float.')

        if not (community_size > 0):
            raise ValueError('`community_size` is expected to be strictly positive.')
        
        if not (len(community_network_degree) == community_size):
            raise ValueError('`community_network_degree` is incompatible with `community_size`.')
        
        if not (all([_ > 0 for _ in community_network_degree])):
            raise ValueError('`community_network_degree` is expected to be strictly positive.')

        if not (0 <= initial_prevalence <= 1):
            raise ValueError('`initial_prevalence` is expected to be of between 0 and 1.')

        if not (0 <= fractional_exponent <= 1):
            raise ValueError('`fractional_exponent` is expected to be of between 0 and 1.')

        for _ in ['mobility rate', 'community infection rate', 'transport infection rate', 'recovery rate', 'immunity loss rate']:
            if not ((_ in epidemic_parameters) and epidemic_parameters[_] >= 0):
                raise ValueError('`epidemic_parameters` is expected to have a field `{}` with non-negative value.'.format(_))


        # SET TRANSPORT NETWORKS

        if not isinstance(transport_networks, (list, tuple)):
            transport_networks = [transport_networks] * 2

        self.transport_networks = transport_networks

        self.fractional_exponent = fractional_exponent

        self.transport_network_adjacency_matrices = [nx.adjacency_matrix(_).toarray().astype('float') for _ in self.transport_networks]

        #def fractional_transition_probability_matrix(A, a):
        #
        #    laplacian = np.diagflat(np.sum(A, axis=1)) - A
        #
        #    if not np.allclose(laplacian, laplacian.T):
        #        raise NotImplementedError('Fractional dynamics cannot be computed for directed networks.')
        #
        #    e, U = np.linalg.eigh(laplacian)
        #    e = np.clip(e, 0, np.trace(laplacian))
        #    e = e**a
        #
        #    laplacian = U @ np.diagflat(e) @ U.T
        #    
        #    transition_matrix = np.eye(*laplacian.shape) -  laplacian / np.diagonal(laplacian)[:,np.newaxis]
        #
        #    return transition_matrix
        #
        #if (self.fractional_exponent != 1):
        #    self.transport_network_transition_matrices = [fractional_transition_probability_matrix(_, self.fractional_exponent) for _ in self.transport_network_adjacency_matrices]
        #else:
        #    self.transport_network_transition_matrices = [_ / np.sum(_, axis=1, keepdims=True) for _ in self.transport_network_adjacency_matrices]
        #self.transport_network_laplacians = [np.eye(*_.shape) - _ for _ in self.transport_network_transition_matrices]


        # SET REMAINING PARAMETERS

        self.world_size = nx.number_of_nodes(self.transport_networks[0])
        if not self.world_size == nx.number_of_nodes(self.transport_networks[1]):
            assert False


        self.transport_network_interpolation_functional = transport_network_interpolation_functional
        
        #self.transport_network_laplacian = lambda t : (1-self.transport_network_interpolation_functional(t)) * self.transport_network_laplacians[0] + self.transport_network_interpolation_functional(t) * self.transport_network_laplacians[1]


        self.community_size = community_size

        self.community_network_degree = community_network_degree
        np.random.shuffle(self.community_network_degree)
        self.community_network_degrees, self.community_network_degree_counts = np.unique(self.community_network_degree, return_counts=True)


        self.initial_site = initial_site
        self.initial_prevalence = initial_prevalence

        self.mobility_rate = epidemic_parameters['mobility rate']
        self.community_infection_rate = epidemic_parameters['community infection rate']
        self.transport_infection_rate = epidemic_parameters['transport infection rate']
        self.recovery_rate = epidemic_parameters['recovery rate']
        self.immunity_loss_rate = epidemic_parameters['immunity loss rate']


        self.epidemic_flavour = self.determine_epidemic_flavour()


        # CREATE A RANDOM EPIDEMIC COMMUNITY NETWORK

        if not skip_community_network_generation:
            print('Generate degree-sequence graph... ', end='')
            self.community_network = nx.random_degree_sequence_graph(self.community_network_degree, tries=100)
            print('Done.')
            self.community_network_adjacency_matrix = nx.adjacency_matrix(self.community_network).toarray().astype('float')
            self.community_contacts = [set(self.community_network[_]) for _ in range(self.community_size)]
        else:
            self.community_network = None
            self.community_network_adjacency_matrix = None
            self.community_contacts = None


        # PRE-INITIALISE STATE

        self.state = {
            'health': np.array(['S' for _ in range(self.community_size)]),
            'site': np.array([0 for _ in range(self.community_size)]),
        }

        self.site_occupancy = [[] for _ in range(self.world_size)]
        self.site_occupancy[0] = list(range(self.community_size))

        self.future_event_times = {
            'move': -np.inf * np.ones(self.community_size),
            'recovery': -np.inf * np.ones(self.community_size),
            'immunity loss': -np.inf * np.ones(self.community_size),
        }


        # CONTROL VARIABLES

        self.event_queue = []

        self.timeseries = {}
        self.simulation_time = np.inf


        self.track_state = False
        self.track_state_motif_counts = False

        #
        self.__reinit__()

    def __reinit__(self):


        def fractional_transition_matrix(A, a):

            laplacian = np.diagflat(np.sum(A, axis=1)) - A

            if (a != 1):
                if not np.allclose(laplacian, laplacian.T):
                    raise NotImplementedError('Fractional dynamics cannot be computed for directed networks.')

                e, U = np.linalg.eigh(laplacian)
                e = np.clip(e, 0, np.trace(laplacian))
                e[e>0] = e[e>0]**a

                laplacian = U @ np.diagflat(e) @ U.T
            
            transition_matrix = np.eye(*laplacian.shape) -  laplacian / np.diagonal(laplacian)[:,np.newaxis]

            return transition_matrix, np.diagonal(laplacian) / np.trace(laplacian)

        self.transport_network_transition_matrices = [None] * len(self.transport_network_adjacency_matrices)
        self.transport_network_equilibrium_distributions = [None] * len(self.transport_network_adjacency_matrices)
        for i, _ in enumerate(self.transport_network_adjacency_matrices):
            self.transport_network_transition_matrices[i], self.transport_network_equilibrium_distributions[i] = fractional_transition_matrix(_, self.fractional_exponent)
        
        self.transport_network_laplacians = [np.eye(*_.shape) - _ for _ in self.transport_network_transition_matrices]


        self.world_size = nx.number_of_nodes(self.transport_networks[0])
        if not self.world_size == nx.number_of_nodes(self.transport_networks[1]):
            assert False

        
        self.transport_network_laplacian = lambda t : (1-self.transport_network_interpolation_functional(t)) * self.transport_network_laplacians[0] + self.transport_network_interpolation_functional(t) * self.transport_network_laplacians[1]


        self.epidemic_flavour = self.determine_epidemic_flavour()


    def __repr__(self):
        return 'Multiplex Transport-Epidemic Simulation ({})\n'.format(self.epidemic_flavour) +\
               '  World size               : {}\n'.format(self.world_size) +\
               '  Community size           : {}\n'.format(self.community_size) +\
               '  Community degree         : {}\n'.format(','.join(['{}*{}'.format(m, k) for k, m in zip(*np.unique(self.community_network_degree, return_counts=True))])) +\
               '  Mobility rate            : {}\n'.format(self.mobility_rate) +\
               '  Community infection rate : {}\n'.format(self.community_infection_rate) +\
               '  Transport infection rate : {}\n'.format(self.transport_infection_rate) +\
               '  Recovery rate            : {}\n'.format(self.recovery_rate) +\
               '  Immunity loss rate       : {}\n'.format(self.immunity_loss_rate) +\
               ''
               
    def initialise_dynamics(self):

        self.event_queue = []
        self.timeseries = {}


        if (0 <= self.initial_site < self.world_size):
            initial_site_distribution = lambda : self.initial_site
        else:
            initial_site_distribution = lambda : rd.randint(0, self.world_size-1)

        for n in range(self.community_size):

            for _ in self.future_event_times.keys():
                self.future_event_times[_][n] = -np.inf

            self.state['health'][n] = 'S'

            x = self.state['site'][n]
            self.site_occupancy[x].remove(n)
            x = initial_site_distribution()
            self.state['site'][n] = x
            self.site_occupancy[x].append(n)


            T = 0.

            if rd.random() < self.initial_prevalence:
                self.event(time=T, subject=n, sender=None, action=self.event.ACTION.COMMUNITY_INFECTION).submit(self.event_queue, lambda e: e.time < self.simulation_time)

            T = 0. + rd.expovariate(1) / self.mobility_rate
            self.future_event_times['move'][n] = T
            self.event(time=T, subject=n, sender=None, action=self.event.ACTION.MOVE).submit(self.event_queue, lambda e: e.time < self.simulation_time)
        

        T = sys.float_info.epsilon
        self.event(time=T, subject=None, sender=None, action=self.event.ACTION.UNSPECIFIED).submit(self.event_queue, lambda e: e.time < self.simulation_time)


    def determine_epidemic_flavour(self):

        def is_zero(x):
            return x == 0

        def is_positive_finite(x):
            return (0 < x and x < np.inf)

        def is_infinite(x):
            return math.isinf(x)

        if is_positive_finite(self.recovery_rate) and is_positive_finite(self.immunity_loss_rate):
            return 'SIRS'
        elif is_zero(self.recovery_rate):
            return 'SI'
        elif is_positive_finite(self.recovery_rate) and is_zero(self.immunity_loss_rate):
            return 'SIR'
        elif is_positive_finite(self.recovery_rate) and is_infinite(self.immunity_loss_rate):
            return 'SIS'
        else:
            warnings.warn('The epidemic parameters do not define a known flavour.', UserWarning)
            return None

    def move_event_handler(self, event):

        r = int(not rd.random() > self.transport_network_interpolation_functional(event.time))


        x = self.state['site'][event.subject]
        self.site_occupancy[x].remove(event.subject)
        x = np.argmax(rd.random() < np.cumsum(self.transport_network_transition_matrices[r][x,:]))
        self.state['site'][event.subject] = x
        self.site_occupancy[x].append(event.subject)


        T = event.time + rd.expovariate(1) / self.mobility_rate
        self.future_event_times['move'][event.subject] = T
        self.event(time=T, subject=event.subject, sender=None, action=self.event.ACTION.MOVE).submit(self.event_queue, lambda e: e.time < self.simulation_time)


        for a in self.site_occupancy[x]:
            if not event.subject == a:
                if self.state['health'][event.subject] == 'I':
                    self.infection_spreading(self.event(time=event.time, subject=a, sender=event.subject, action=self.event.ACTION.TRANSPORT_INFECTION))

                if self.state['health'][a] == 'I':
                    self.infection_spreading(self.event(time=event.time, subject=event.subject, sender=a, action=self.event.ACTION.TRANSPORT_INFECTION))

    def infection_event_handler(self, event):

        if (self.state['health'][event.subject] == 'S'):
            
            self.state['health'][event.subject] = 'I'
            
            
            # DETERMINE FUTURE RECOVERY (I --> R) TIME
            T = event.time + rd.expovariate(1) * (1/self.recovery_rate if self.recovery_rate != 0 else np.inf)
            self.future_event_times['recovery'][event.subject] = T
            self.event(time=T, subject=event.subject, sender=None, action=self.event.ACTION.RECOVERY).submit(self.event_queue, lambda e: e.time < self.simulation_time)
            
            # DETERMINE FUTURE IMMUNITY LOSS (R --> S) TIME
            T = T + rd.expovariate(1) * (1/self.immunity_loss_rate if self.immunity_loss_rate != 0 else np.inf)
            self.future_event_times['immunity loss'][event.subject] = T
            self.event(time=T, subject=event.subject, sender=None, action=self.event.ACTION.IMMUNITY_LOSS).submit(self.event_queue, lambda e: e.time < self.simulation_time)

            
            for n in self.community_contacts[event.subject]:
                self.infection_spreading(self.event(time=event.time, subject=n, sender=event.subject, action=self.event.ACTION.COMMUNITY_INFECTION))

            for n in self.site_occupancy[self.state['site'][event.subject]]:
                if not event.subject == n:
                    self.infection_spreading(self.event(time=event.time, subject=n, sender=event.subject, action=self.event.ACTION.TRANSPORT_INFECTION))

        else:
            self.infection_spreading(event)

    def infection_spreading(self, event):

        if event.action == '':
            pass

        elif event.action == self.event.ACTION.COMMUNITY_INFECTION:
            infection_rate = self.community_infection_rate
            T_infection_path_breaking = np.inf

        elif event.action == self.event.ACTION.TRANSPORT_INFECTION:
            infection_rate = self.transport_infection_rate
            T_infection_path_breaking = min(self.future_event_times['move'][event.subject], self.future_event_times['move'][event.sender])

        else:
            pass


        if (event.action == self.event.ACTION.COMMUNITY_INFECTION or event.action == self.event.ACTION.TRANSPORT_INFECTION):
            if self.future_event_times['immunity loss'][event.subject] < self.future_event_times['recovery'][event.sender]:
                T = max(event.time, self.future_event_times['immunity loss'][event.subject]) + rd.expovariate(1) * (1/infection_rate if infection_rate != 0 else np.inf)

                if (T < self.future_event_times['recovery'][event.sender] and T < T_infection_path_breaking):
                    if T < self.simulation_time:
                        self.event(time=T, subject=event.subject, sender=event.sender, action=event.action).submit(self.event_queue, lambda e: e.time < self.simulation_time)

    def recovery_event_handler(self, event):

        self.state['health'][event.subject] = 'R'
        self.future_event_times['recovery'][event.subject] = -np.inf

    def immunity_loss_event_handler(self, event):

        self.state['health'][event.subject] = 'S'
        self.future_event_times['immunity loss'][event.subject] = -np.inf


    def update_timeseries(self, time):

        __observable_values = {}

        for h in {'S', 'I', 'R'}:
            __observable_values[h] = np.mean(self.state['health'] == h)

        if self.track_state_motif_counts:
            for motif, count in self.count_state_motifs().items():
                __observable_values['#{}'.format(motif)] = count

        if self.track_state:
            __observable_values['state'] = deepcopy(self.state)

        
        #self.timeseries = self.timeseries.append(
        #    pd.Series(__observable_values, name=time)
        #)

        if not self.timeseries:
            self.timeseries['t'] = [time]
            for observable, value in __observable_values.items():
                self.timeseries[observable] = [value]
        else:
            self.timeseries['t'].append(time)
            for observable, value in __observable_values.items():
                self.timeseries[observable].append(value)


    def count_state_motifs(self, include_triples=False):

        A_ = {
            'c': self.community_network_adjacency_matrix.astype('int'),
            't': (self.state['site'][np.newaxis,:] == self.state['site'][:,np.newaxis]) * (1 - np.eye(self.community_size, dtype='int')),
        }

        for w in A_.keys():
            A_[w] = scipy.sparse.csr_array(A_[w])
        
        E_ = {
            _: scipy.sparse.dia_array(((self.state['health'] == _).astype('int'), 0), shape=(self.community_size,self.community_size)) for _ in ['S', 'I', 'R']
        }

        def H(h):
            C = E_[h]
            return (C.sum())

        def HwH(hh, w):
            #return np.dot((H == hh[1]), A[w] @ (H == hh[0]))
            C = E_[hh[0]] @ A_[w] @ E_[hh[1]]
            return (C.sum())

        def HwHwH(hhh, ww):
            C = E_[hhh[0]] @ A_[ww[0]] @ E_[hhh[1]] @ A_[ww[1]] @ E_[hhh[2]]
            return (C.sum() - (ww[0] == ww[1]) * C.trace())
        

        S = H('S')
        I = H('I')
        R = H('R')


        ScS = HwH(('S', 'S'), 'c')
        IcI = HwH(('I', 'I'), 'c')
        RcR = HwH(('R', 'R'), 'c')

        ScI = HwH(('S', 'I'), 'c')
        ScR = HwH(('S', 'R'), 'c')
        IcR = HwH(('I', 'R'), 'c')

        StS = HwH(('S', 'S'), 't')
        ItI = HwH(('I', 'I'), 't')
        RtR = HwH(('R', 'R'), 't')

        StI = HwH(('S', 'I'), 't')
        StR = HwH(('S', 'R'), 't')
        ItR = HwH(('I', 'R'), 't')


        if include_triples:
            ScScI = HwHwH(('S', 'S', 'I'), ('c', 'c'))
            ScStI = HwHwH(('S', 'S', 'I'), ('c', 't'))
            StScI = HwHwH(('S', 'S', 'I'), ('t', 'c'))
            StStI = HwHwH(('S', 'S', 'I'), ('t', 't'))
            IcScI = HwHwH(('I', 'S', 'I'), ('c', 'c'))
            IcStI = HwHwH(('I', 'S', 'I'), ('c', 't'))
            ItStI = HwHwH(('I', 'S', 'I'), ('t', 't'))
            RcScI = HwHwH(('R', 'S', 'I'), ('c', 'c'))
            RcStI = HwHwH(('R', 'S', 'I'), ('c', 't'))
            RtScI = HwHwH(('R', 'S', 'I'), ('t', 'c'))
            RtStI = HwHwH(('R', 'S', 'I'), ('t', 't'))


        if include_triples:
            motifs = (S, I, R, ScS, IcI, RcR, ScI, ScR, IcR, StS, ItI, RtR, StI, StR, ItR, ScScI, ScStI, StScI, StStI, IcScI, IcStI, ItStI, RcScI, RcStI, RtScI, RtStI)
            motif_labels = ('S', 'I', 'R', 'ScS', 'IcI', 'RcR', 'ScI', 'ScR', 'IcR', 'StS', 'ItI', 'RtR', 'StI', 'StR', 'ItR', 'ScScI', 'ScStI', 'StScI', 'StStI', 'IcScI', 'IcStI', 'ItStI', 'RcScI', 'RcStI', 'RtScI', 'RtStI')
        else:
            motifs = (S, I, R, ScS, IcI, RcR, ScI, ScR, IcR, StS, ItI, RtR, StI, StR, ItR)
            motif_labels = ('S', 'I', 'R', 'ScS', 'IcI', 'RcR', 'ScI', 'ScR', 'IcR', 'StS', 'ItI', 'RtR', 'StI', 'StR', 'ItR')

        if False:
            l = max(math.ceil(math.log10(max(*motifs))) + 1, 6)

            print((('{: >' + '{:d}'.format(l) + 's}') * len(motifs)).format(*motif_labels))
            print((('{: >' + '{:d}'.format(l) + 'd}') * len(motifs)).format(*motifs))


        return dict(zip(motif_labels, motifs))

    #def resampled_timeseries(self, T):
    #
    #    __timeseries = self.timeseries
    #
    #    __NaN_series = pd.DataFrame(columns=__timeseries.columns, index=T)
    #    interpolated_timeseries = __timeseries.combine_first(__NaN_series).fillna(method='ffill')
    #
    #    return pd.concat([interpolated_timeseries[c][T] for c in __timeseries.columns], axis=1)


    def simulate(self, time, track_state=False, track_state_motif_counts=False):

        if self.community_network is None:
            raise RuntimeError('Simulation is only possible with a community network.')

        self.track_state = track_state
        self.track_state_motif_counts = track_state_motif_counts

        self.simulation_time = time
        self.initialise_dynamics()


        while len(self.event_queue) > 0:

            event = heapq.heappop(self.event_queue)
            #print(event.time, event.subject, event.action)
            if event.time > self.simulation_time:
                break

            if event.action == self.event.ACTION.MOVE:
                self.move_event_handler(event)
            elif event.action == self.event.ACTION.COMMUNITY_INFECTION or event.action == self.event.ACTION.TRANSPORT_INFECTION:
                self.infection_event_handler(event)
            elif event.action == self.event.ACTION.RECOVERY:
                self.recovery_event_handler(event)
            elif event.action == self.event.ACTION.IMMUNITY_LOSS:
                self.immunity_loss_event_handler(event)
            else:
                pass

            if event.time > 0 and event.action != self.event.ACTION.MOVE:
                self.update_timeseries(event.time)


        self.update_timeseries(self.simulation_time)

        return pd.DataFrame.from_dict(self.timeseries).set_index('t')


    def macroscopic_mean_field_timeseries(self, T, order=1, reduced=False, normalised_output=True, full_output=True):

        self.epidemic_flavour = self.determine_epidemic_flavour()



        if len(self.community_network_degrees) > 1:

            if self.epidemic_flavour in ['SIRS', 'SI', 'SIR']:
                return self.SIRS_irregular_macroscopic_mean_field_timeseries(T, order=order, normalised_output=normalised_output, full_output=full_output)

            if self.epidemic_flavour == 'SIS':
                return self.SIS_irregular_macroscopic_mean_field_timeseries(T, order=order, normalised_output=normalised_output, full_output=full_output)

        else:

            ## THIS IS A MEAN FIELD MODEL VALID ONLY FOR k-REGULAR (EPIDEMIC) COMMUNITY NETWORKS
            #
            #k = np.mean([v for _, v in self.community_network.degree()])
            #if np.std([v for _, v in self.community_network.degree()]) > 0:
            #    warnings.warn('The implemented macroscopic mean field model is not valid for non-regular community networks.', UserWarning)
            #
            #assert all([_ == k for _ in self.community_network_degree])


            if self.epidemic_flavour in ['SIRS', 'SI', 'SIR']:
                return self.SIRS_regular_macroscopic_mean_field_timeseries(T, order=order, reduced=reduced, normalised_output=normalised_output, full_output=full_output)

            if self.epidemic_flavour == 'SIS':
                return self.SIS_regular_macroscopic_mean_field_timeseries(T, order=order, reduced=reduced, normalised_output=normalised_output, full_output=full_output)


        raise NotImplementedError('An {} mean-field model has not been implemented yet.'.format(self.epidemic_flavour))


    def SIRS_regular_macroscopic_mean_field_timeseries(self, T, order=1, reduced=True, normalised_output=True, full_output=True):
        
        tolerance = 1e-10 # used to determine whether initial conditions are compatible with the reduced model

        # INITIAL CONDITIONS

        if (0 <= self.initial_site < self.world_size):
            p0 = np.zeros(self.world_size)
            p0[self.initial_site] = 1
        else:
            p0 = np.ones(self.world_size) / self.world_size


        c0 = np.dot(p0, p0)
        S0 = (1 - self.initial_prevalence) * self.community_size
        I0 = self.initial_prevalence * self.community_size
        R0 = 0

        ScS0 = self.community_network_degree[0] / self.community_size * S0**2
        IcI0 = self.community_network_degree[0] / self.community_size * I0**2
        RcR0 = self.community_network_degree[0] / self.community_size * R0**2
        ScI0 = self.community_network_degree[0] / self.community_size * S0 * I0
        ScR0 = self.community_network_degree[0] / self.community_size * S0 * R0
        IcR0 = self.community_network_degree[0] / self.community_size * I0 * R0
        
        StS0 = c0 * S0**2
        ItI0 = c0 * I0**2
        RtR0 = c0 * R0**2
        StI0 = c0 * S0 * I0
        StR0 = c0 * S0 * R0
        ItR0 = c0 * I0 * R0


        if order == 1:

            if not reduced:
                
                def dV_dt(t, V, mu, beta_c, beta_t, gamma, sigma, N, k):

                    p = V[0:self.world_size]
                    S, I, R = V[self.world_size:]
    
                    c = np.dot(p, p)
    
                    ScI = (k / N) * I * S
                    StI = c * I * S


                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p

                    dS_dt = -(beta_c * ScI + beta_t * StI) + sigma * R
                    dI_dt = (beta_c * ScI + beta_t * StI) - gamma * I
                    dR_dt = gamma * I - sigma * R
    
                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dI_dt, dR_dt])))
    
                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0, I0, R0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)
    
                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, I_t, R_t = (_ for _ in V_t[self.world_size:,:])

            else:

                if not abs(S0 + I0 + R0 - self.community_size) < tolerance:
                    warnings.warn('Initial conditions are not consistent with the reduced model (Tolerance: {}).'.format(tolerance))

                def dV_dt(t, V, mu, beta_c, beta_t, gamma, sigma, N, k):

                    p = V[0:self.world_size]
                    S, I = V[self.world_size:]
    
                    c = np.dot(p, p)

                    R = N - (S + I)

    
                    ScI = (k / N) * I * S
                    StI = c * I * S


                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p

                    dS_dt = -(beta_c * ScI + beta_t * StI) + sigma * R
                    dI_dt = (beta_c * ScI + beta_t * StI) - gamma * I
    
                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dI_dt])))
    
                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0, I0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)
    
                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, I_t = (_ for _ in V_t[self.world_size:,:])

                R_t = self.community_size - (S_t + I_t)

            if normalised_output:
                S_t /= self.community_size
                I_t /= self.community_size
                R_t /= self.community_size

            if full_output:
                return pd.DataFrame(data=dict(zip(('S', 'I', 'R',), (S_t, I_t, R_t,))), index=t)
            else:
                return pd.DataFrame(data=dict(zip(('I',), (I_t,))), index=t)

        elif order == 2:

            if not reduced:
                
                def dV_dt(t, V, mu, beta_c, beta_t, gamma, sigma, N, k):
                    
                    if gamma == 0:
                        sigma = 0

                    p = V[0:self.world_size]
                    S, I, R, ScS, IcI, RcR, ScI, ScR, IcR, StS, ItI, RtR, StI, StR, ItR = V[self.world_size:]

                    c = np.dot(p, p)

                    if S > 0:
                        ScScI = (1 - 1/k) * ScS * ScI / S
                        ScStI = ScS * StI / S
                        IcScI = (1 - 1/k) * ScI**2 / S
                        IcStI = ScI * StI / S
                        RcScI = (1 - 1/k) * ScR * ScI / S
                        RcStI = ScR * StI / S

                        StScI = StS * ScI / S
                        StStI = (1 - 1/(c * N)) * StS * StI / S
                        ItScI = StI * ScI / S
                        ItStI = (1 - 1/(c * N)) * StI**2 / S
                        RtScI = StR * ScI / S
                        RtStI = (1 - 1/(c * N)) * StR * StI / S
                    else:
                        ScScI = 0
                        ScStI = 0
                        IcScI = 0
                        IcStI = 0
                        RcScI = 0
                        RcStI = 0

                        StScI = 0
                        StStI = 0
                        ItScI = 0
                        ItStI = 0
                        RtScI = 0
                        RtStI = 0

                    
                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                    dc_dt = 2 * np.dot(p, dp_dt)

                    dS_dt = sigma * R - beta_c * ScI - beta_t * StI
                    dI_dt = beta_c * ScI + beta_t * StI - gamma * I
                    dR_dt = gamma * I - sigma * R

                    dScS_dt = -2 * (beta_c * ScScI + beta_t * ScStI) + 2 * sigma * ScR
                    dIcI_dt = 2 * (beta_c * IcScI + beta_t * IcStI) + 2 * beta_c * ScI - 2 * gamma * IcI
                    dRcR_dt = 2 * gamma * IcR - 2 * sigma * RcR
                    dScI_dt = (beta_c * (ScScI - IcScI) + beta_t * (ScStI - IcStI)) - (beta_c + gamma) * ScI + sigma * IcR
                    dScR_dt = - (beta_c * RcScI + beta_t * RcStI) + gamma * ScI + sigma * (RcR - ScR)
                    dIcR_dt = (beta_c * RcScI + beta_t * RcStI) + gamma * IcI - (gamma + sigma) * IcR

                    dStS_dt = -2 * (beta_c * StScI + beta_t * StStI) + 2 * sigma * StR \
                              + dc_dt * S**2
                    dItI_dt = 2 * (beta_c * ItScI + beta_t * ItStI) + 2 * beta_t * StI - 2 * gamma * ItI \
                              + dc_dt * I**2
                    dRtR_dt = 2 * gamma * ItR - 2 * sigma * RtR \
                              + dc_dt * R**2
                    dStI_dt = (beta_c * (StScI - ItScI) + beta_t * (StStI - ItStI)) - (beta_t + gamma) * StI + sigma * ItR \
                              + dc_dt * S * I
                    dStR_dt = - (beta_c * RtScI + beta_t * RtStI) + gamma * StI + sigma * (RtR - StR) \
                              + dc_dt * S * R
                    dItR_dt = (beta_c * RtScI + beta_t * RtStI) + gamma * ItI - (gamma + sigma) * ItR \
                              + dc_dt * I * R

                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dI_dt, dR_dt, dScS_dt, dIcI_dt, dRcR_dt, dScI_dt, dScR_dt, dIcR_dt, dStS_dt, dItI_dt, dRtR_dt, dStI_dt, dStR_dt, dItR_dt])))

                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0, I0, R0, ScS0, IcI0, RcR0, ScI0, ScR0, IcR0, StS0, ItI0, RtR0, StI0, StR0, ItR0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, I_t, R_t, ScS_t, IcI_t, RcR_t, ScI_t, ScR_t, IcR_t, StS_t, ItI_t, RtR_t, StI_t, StR_t, ItR_t = (_ for _ in V_t[self.world_size:,:])

            else:

                if not(abs(S0 + I0 + R0 - self.community_size) < tolerance and \
                       abs(2 * (ScI0 + IcR0 + ScR0) + ScS0 + IcI0 + RcR0 - self.community_network_degree[0] * self.community_size) < tolerance and \
                       abs(2 * (StI0 + ItR0 + StR0) + StS0 + ItI0 + RtR0 - c0 * self.community_size**2) < tolerance and \
                       abs(ScS0 + ScI0 + ScR0 - self.community_network_degree[0] * S0) < tolerance and \
                       abs(ScI0 + IcI0 + IcR0 - self.community_network_degree[0] * I0) < tolerance and \
                       abs(ScR0 + IcR0 + RcR0 - self.community_network_degree[0] * I0) < tolerance and \
                       abs(StS0 + StI0 + StR0 - c0 * self.community_size * S0) < tolerance and \
                       abs(StI0 + ItI0 + ItR0 - c0 * self.community_size * I0) < tolerance and \
                       abs(StR0 + ItR0 + RtR0 - c0 * self.community_size * I0) < tolerance \
                ):
                    warnings.warn('Initial conditions are not consistent with the reduced model (Tolerance: {}).'.format(tolerance))
                
                def dV_dt(t, V, mu, beta_c, beta_t, gamma, sigma, N, k):
                    
                    if gamma == 0:
                        sigma = 0

                    p = V[0:self.world_size]
                    S, I, ScS, IcI, ScI, StS, ItI, StI = V[self.world_size:]

                    c = np.dot(p, p)

                    R = N - (S + I)
                    ScR = k * S - ScS - ScI
                    IcR = k * I - IcI - ScI
                    RcR = k * R - ScR - IcR
                    StR = c * N * S - StS - StI
                    ItR = c * N * I - ItI - StI
                    RtR = c * N * R - StR - ItR


                    if S > 0:
                        ScScI = (1 - 1/k) * ScS * ScI / S
                        ScStI = ScS * StI / S
                        IcScI = (1 - 1/k) * ScI**2 / S
                        IcStI = ScI * StI / S
                        RcScI = (1 - 1/k) * ScR * ScI / S
                        RcStI = ScR * StI / S

                        StScI = StS * ScI / S
                        StStI = (1 - 1/(c * N)) * StS * StI / S
                        ItScI = StI * ScI / S
                        ItStI = (1 - 1/(c * N)) * StI**2 / S
                        RtScI = StR * ScI / S
                        RtStI = (1 - 1/(c * N)) * StR * StI / S
                    else:
                        ScScI = 0
                        ScStI = 0
                        IcScI = 0
                        IcStI = 0
                        RcScI = 0
                        RcStI = 0

                        StScI = 0
                        StStI = 0
                        ItScI = 0
                        ItStI = 0
                        RtScI = 0
                        RtStI = 0

                    
                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                    dc_dt = 2 * np.dot(p, dp_dt)

                    dS_dt = sigma * R - beta_c * ScI - beta_t * StI
                    dI_dt = beta_c * ScI + beta_t * StI - gamma * I

                    dScS_dt = -2 * (beta_c * ScScI + beta_t * ScStI) + 2 * sigma * ScR
                    dIcI_dt = 2 * (beta_c * IcScI + beta_t * IcStI) + 2 * beta_c * ScI - 2 * gamma * IcI
                    dScI_dt = (beta_c * (ScScI - IcScI) + beta_t * (ScStI - IcStI)) - (beta_c + gamma) * ScI + sigma * IcR

                    dStS_dt = -2 * (beta_c * StScI + beta_t * StStI) + 2 * sigma * StR \
                              + dc_dt * S**2
                    dItI_dt = 2 * (beta_c * ItScI + beta_t * ItStI) + 2 * beta_t * StI - 2 * gamma * ItI \
                              + dc_dt * I**2
                    dStI_dt = (beta_c * (StScI - ItScI) + beta_t * (StStI - ItStI)) - (beta_t + gamma) * StI + sigma * ItR \
                              + dc_dt * S * I

                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dI_dt, dScS_dt, dIcI_dt, dScI_dt, dStS_dt, dItI_dt, dStI_dt])))

                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0, I0, ScS0, IcI0, ScI0, StS0, ItI0, StI0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, I_t, ScS_t, IcI_t, ScI_t, StS_t, ItI_t, StI_t = (_ for _ in V_t[self.world_size:,:])

                R_t = self.community_size - (S_t + I_t)
                ScR_t = self.community_network_degree[0] * S_t - ScS_t - ScI_t
                IcR_t = self.community_network_degree[0] * I_t - IcI_t - ScI_t
                RcR_t = self.community_network_degree[0] * R_t - ScR_t - IcR_t
                StR_t = c_t * self.community_size * S_t - StS_t - StI_t
                ItR_t = c_t * self.community_size * I_t - ItI_t - StI_t
                RtR_t = c_t * self.community_size * R_t - StR_t - ItR_t

            if normalised_output:
                S_t /= self.community_size
                I_t /= self.community_size
                R_t /= self.community_size

                ScS_t /= self.community_size * (self.community_size - 1)
                IcI_t /= self.community_size * (self.community_size - 1)
                RcR_t /= self.community_size * (self.community_size - 1)
                ScI_t /= self.community_size * (self.community_size - 1)
                ScR_t /= self.community_size * (self.community_size - 1)
                IcR_t /= self.community_size * (self.community_size - 1)
                StS_t /= self.community_size * (self.community_size - 1)
                ItI_t /= self.community_size * (self.community_size - 1)
                RtR_t /= self.community_size * (self.community_size - 1)
                StI_t /= self.community_size * (self.community_size - 1)
                StR_t /= self.community_size * (self.community_size - 1)
                ItR_t /= self.community_size * (self.community_size - 1)

            if full_output:
                return pd.DataFrame(data=dict(zip(('S', 'I', 'R', 'ScS', 'IcI', 'RcR', 'ScI', 'ScR', 'IcR', 'StS', 'ItI', 'RtR', 'StI', 'StR', 'ItR',), (S_t, I_t, R_t, ScS_t, IcI_t, RcR_t, ScI_t, ScR_t, IcR_t, StS_t, ItI_t, RtR_t, StI_t, StR_t, ItR_t,))), index=t)
            else:
                return pd.DataFrame(data=dict(zip(('I',), (I_t,))), index=t)


    def SIS_regular_macroscopic_mean_field_timeseries(self, T, order=1, reduced=True, normalised_output=True, full_output=True):

        tolerance = 1e-10 # used to determine whether initial conditions are compatible with the reduced model

        # INITIAL CONDITIONS

        if (0 <= self.initial_site < self.world_size):
            p0 = np.zeros(self.world_size)
            p0[self.initial_site] = 1
        else:
            p0 = np.ones(self.world_size) / self.world_size


        c0 = np.dot(p0, p0)
        S0 = (1 - self.initial_prevalence) * self.community_size
        I0 = self.initial_prevalence * self.community_size

        ScI0 = self.community_network_degree[0] / self.community_size * S0 * I0
        ScS0 = self.community_network_degree[0] / self.community_size * S0**2
        IcI0 = self.community_network_degree[0] / self.community_size * I0**2
        
        StI0 = c0 * S0 * I0
        StS0 = c0 * S0**2
        ItI0 = c0 * I0**2


        if order == 1:

            if not reduced:

                def dV_dt(t, V, mu, beta_c, beta_t, gamma, k):
                
                    p = V[0:self.world_size]
                    S, I = V[self.world_size:]

                    N = self.community_size
    
                    c = np.dot(p, p)
    
                    ScI = (k / N) * I * S
                    StI = c * I * S


                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p

                    dS_dt = -(beta_c * ScI + beta_t * StI) + gamma * I
                    dI_dt = (beta_c * ScI + beta_t * StI) - gamma * I
    
                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dI_dt])))
    
                    return __dV_dt
    
                V0 = np.concatenate((p0, np.array([S0, I0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)
    
                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, I_t = (_ for _ in V_t[self.world_size:,:])

            else:

                if not abs(S0 + I0 - self.community_size) < tolerance:
                    warnings.warn('Initial conditions are not consistent with the reduced model (Tolerance: {}).'.format(tolerance))

                def dV_dt(t, V, mu, beta_c, beta_t, gamma, N, k):

                    p = V[0:self.world_size]
                    S = V[self.world_size]

                    c = np.dot(p, p)

                    I = N - S


                    ScI = (k / N) * I * S
                    StI = c * I * S

                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                    dS_dt = -(beta_c * ScI + beta_t * StI) + gamma * I

                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt])))

                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t = V_t[self.world_size,:]

                I_t = self.community_size - S_t

            if normalised_output:
                S_t /= self.community_size
                I_t /= self.community_size

            if full_output:
                #return t, S_t, I_t
                return pd.DataFrame(data=dict(zip(('S', 'I',), (S_t, I_t,))), index=t)
            else:
                #return t, I_t
                return pd.DataFrame(data=dict(zip(('I',), (I_t,))), index=t)

        elif order == 2:

            if not reduced:

                def dV_dt(t, V, mu, beta_c, beta_t, gamma, N, k):

                    p = V[0:self.world_size]
                    S, I, ScI, ScS, IcI, StI, StS, ItI = V[self.world_size:]

                    c = np.dot(p, p)

                    if S > 0:
                        ScScI = (1 - 1/k) * ScS * ScI / S
                        StScI = StS * ScI / S
                        ScStI = ScS * StI / S
                        StStI = (1 - 1/(c * N)) * StS * StI / S
                        IcScI = (1 - 1/k) * ScI**2 / S
                        IcStI = ScI * StI / S
                        ItStI = (1 - 1/(c * N)) * StI**2 / S
                    else:
                        ScScI = 0
                        StScI = 0
                        ScStI = 0
                        StStI = 0
                        IcScI = 0
                        IcStI = 0
                        ItStI = 0


                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                    dc_dt = 2 * np.dot(p, dp_dt)

                    dS_dt = gamma * I - beta_c * ScI - beta_t * StI
                    dI_dt = beta_c * ScI + beta_t * StI - gamma * I

                    dScI_dt = (beta_c * (ScScI - IcScI) + beta_t * (ScStI - IcStI)) - beta_c * ScI - gamma * (ScI - IcI)
                    dScS_dt = -2 * (beta_c * ScScI + beta_t * ScStI) + 2 * gamma * ScI
                    dIcI_dt = 2 * (beta_c * IcScI + beta_t * IcStI) + 2 * beta_c * ScI - 2 * gamma * IcI

                    dStI_dt = (beta_c * (StScI - IcStI) + beta_t * (StStI - ItStI)) - beta_t * StI - gamma * (StI - ItI) \
                              + dc_dt * S * I
                    dStS_dt = -2 * (beta_c * StScI + beta_t * StStI) + 2 * gamma * StI \
                              + dc_dt * S**2
                    dItI_dt = 2 * (beta_c * IcStI + beta_t * ItStI) + 2 * beta_t * StI - 2 * gamma * ItI \
                              + dc_dt * I**2

                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dI_dt, dScI_dt, dScS_dt, dIcI_dt, dStI_dt, dStS_dt, dItI_dt])))

                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0, I0, ScI0, ScS0, IcI0, StI0, StS0, ItI0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, I_t, ScI_t, ScS_t, IcI_t, StI_t, StS_t, ItI_t = (_ for _ in V_t[self.world_size:,:])
                
            else:

                if not(abs(S0 + I0 - self.community_size) < tolerance and \
                       abs(2 * ScI0 + ScS0 + IcI0 - self.community_network_degree[0] * self.community_size) < tolerance and \
                       abs(2 * StI0 + StS0 + ItI0 - c0 * self.community_size**2) < tolerance and \
                       abs(ScS0 + ScI0 - self.community_network_degree[0] * S0) < tolerance and \
                       abs(ScI0 + IcI0 - self.community_network_degree[0] * I0) < tolerance and \
                       abs(StS0 + StI0 - c0 * self.community_size * S0) < tolerance and \
                       abs(StI0 + ItI0 - c0 * self.community_size * I0) < tolerance \
                ):
                    warnings.warn('Initial conditions are not consistent with the reduced model (Tolerance: {}).'.format(tolerance))

                def dV_dt(t, V, mu, beta_c, beta_t, gamma, N, k):

                    p = V[0:self.world_size]
                    S, ScS, StS = V[self.world_size:]

                    c = np.dot(p, p)

                    I = N - S
                    ScI = k * S - ScS
                    IcI = k * I - ScI
                    StI = c * N * S - StS
                    ItI = c * N * I - StI


                    if S > 0:
                        ScScI = (1 - 1/k) * ScS * ScI / S
                        StScI = StS * ScI / S
                        ScStI = ScS * StI / S
                        StStI = (1 - 1/(c * N)) * StS * StI / S
                        IcScI = (1 - 1/k) * ScI**2 / S
                        IcStI = ScI * StI / S
                        ItStI = (1 - 1/(c * N)) * StI**2 / S
                    else:
                        ScScI = 0
                        StScI = 0
                        ScStI = 0
                        StStI = 0
                        IcScI = 0
                        IcStI = 0
                        ItStI = 0
                    

                    dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                    dc_dt = 2 * np.dot(p, dp_dt)

                    dS_dt = gamma * I - beta_c * ScI - beta_t * StI

                    dScS_dt = -2 * (beta_c * ScScI + beta_t * ScStI) + 2 * gamma * ScI

                    dStS_dt = -2 * (beta_c * StScI + beta_t * StStI) + 2 * gamma * StI \
                              + dc_dt * S**2

                    __dV_dt = np.concatenate((dp_dt, np.array([dS_dt, dScS_dt, dStS_dt])))

                    return __dV_dt

                V0 = np.concatenate((p0, np.array([S0, ScS0, StS0])))
                solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_size, self.community_network_degree[0]), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

                t = solve_ivp_solution.t
                V_t = solve_ivp_solution.y

                c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
                S_t, ScS_t, StS_t = (_ for _ in V_t[self.world_size:,:])

                I_t = self.community_size - S_t
                ScI_t = self.community_network_degree[0] * S_t - ScS_t
                IcI_t = self.community_network_degree[0] * I_t - ScI_t
                StI_t = c_t * self.community_size * S_t - StS_t
                ItI_t = c_t * self.community_size * I_t - StI_t


            if normalised_output:
                S_t /= self.community_size
                I_t /= self.community_size

                ScI_t /= self.community_size * (self.community_size - 1)
                ScS_t /= self.community_size * (self.community_size - 1)
                IcI_t /= self.community_size * (self.community_size - 1)
                StI_t /= self.community_size * (self.community_size - 1)
                StS_t /= self.community_size * (self.community_size - 1)
                ItI_t /= self.community_size * (self.community_size - 1)

            if full_output:
                #return t, S_t, I_t, ScI_t, ScS_t, IcI_t, StI_t, StS_t, ItI_t
                return pd.DataFrame(data=dict(zip(('S', 'I', 'ScI', 'ScS', 'IcI', 'StI', 'StS', 'ItI',), (S_t, I_t, ScI_t, ScS_t, IcI_t, StI_t, StS_t, ItI_t,))), index=t)
            else:
                #return t, I_t
                return pd.DataFrame(data=dict(zip(('I',), (I_t,))), index=t)


    def SIRS_irregular_macroscopic_mean_field_timeseries(self, T, order=1, normalised_output=True, full_output=True):

        # INITIAL CONDITIONS

        if (0 <= self.initial_site < self.world_size):
            p0 = np.zeros(self.world_size)
            p0[self.initial_site] = 1
        else:
            p0 = np.ones(self.world_size) / self.world_size


        c0 = np.dot(p0, p0)
        S_0 = (1 - self.initial_prevalence) * self.community_network_degree_counts
        I_0 = self.initial_prevalence * self.community_network_degree_counts
        R_0 = 0 * self.community_network_degree_counts


        S_cS_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(S_0, S_0)
        I_cI_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(I_0, I_0)
        R_cR_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(R_0, R_0)
        S_cI_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(S_0, I_0)
        S_cR_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(S_0, R_0)
        I_cR_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(I_0, R_0)
        
        S_tS_0 = c0 * np.outer(S_0, S_0)
        I_tI_0 = c0 * np.outer(I_0, I_0)
        R_tR_0 = c0 * np.outer(R_0, R_0)
        S_tI_0 = c0 * np.outer(S_0, I_0)
        S_tR_0 = c0 * np.outer(S_0, R_0)
        I_tR_0 = c0 * np.outer(I_0, R_0)


        if order == 1:

            def dV_dt(t, V, mu, beta_c, beta_t, gamma, sigma, N, K):

                p = V[0:self.world_size]

                M = len(K)
                S_, I_, R_, _ = np.vsplit(V[self.world_size:].reshape(1 * 3,M), np.cumsum([1] * 3))
    
                c = np.dot(p, p)
    
                S_cI = (np.outer(K, K) / np.dot(K, self.community_network_degree_counts) * np.outer(S_, I_)).sum(axis=1)
                S_tI = (c * np.outer(S_, I_)).sum(axis=1)


                dp_dt = - mu * self.transport_network_laplacian(t).T @ p

                dS__dt = -(beta_c * S_cI + beta_t * S_tI) + sigma * R_
                dI__dt = (beta_c * S_cI + beta_t * S_tI) - gamma * I_
                dR__dt = gamma * I_ - sigma * R_
    
                __dV_dt = np.concatenate((dp_dt, np.vstack([dS__dt, dI__dt, dR__dt]).reshape(-1)))
    
                return __dV_dt

            V0 = np.concatenate((p0, np.vstack([S_0, I_0, R_0]).reshape(-1)))
            solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degrees), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)
    
            t = solve_ivp_solution.t
            V_t = solve_ivp_solution.y

            c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
            M = len(self.community_network_degrees)
            S__t, I__t, R__t, _ = np.vsplit(V_t[self.world_size:,:].reshape(1 * 3,M,len(T)), np.cumsum([1] * 3))
            S__t = S__t.squeeze()
            I__t = I__t.squeeze()
            R__t = R__t.squeeze()

            
            if normalised_output:
                S__t /= self.community_size
                I__t /= self.community_size
                R__t /= self.community_size

            if full_output:
                return pd.DataFrame(data=dict(zip(('S', 'I', 'R',), (S__t.sum(axis=0), I__t.sum(axis=0), R__t.sum(axis=0),))), index=t)
            else:
                return pd.DataFrame(data=dict(zip(('I',), (I__t.sum(axis=0),))), index=t)

        elif order == 2:
            
            def dV_dt(t, V, mu, beta_c, beta_t, gamma, sigma, N, K):
                
                if gamma == 0:
                    sigma = 0

                p = V[0:self.world_size]
                
                M = len(K)
                S_, I_, R_, S_cS_, I_cI_, R_cR_, S_cI_, S_cR_, I_cR_, S_tS_, I_tI_, R_tR_, S_tI_, S_tR_, I_tR_, _ = np.vsplit(V[self.world_size:].reshape(1 * 3 + M * 2*6,M), np.cumsum([1] * 3 + [M] * 2*6))
                S_ = S_.flatten()
                I_ = I_.flatten()
                R_ = R_.flatten()


                c = np.dot(p, p)


                S_cI = S_cI_.sum(axis=1)
                S_tI = S_tI_.sum(axis=1)

                S_cS_cI = S_cS_ * np.where(S_ > 0, (1 - 1/K) * S_cI / S_, 0)[np.newaxis,:]
                S_cS_tI = S_cS_ * np.where(S_ > 0, S_tI / S_, 0)[np.newaxis,:]
                I_cS_cI = S_cI_.T * np.where(S_ > 0, (1 - 1/K) * S_cI / S_, 0)[np.newaxis,:]
                I_cS_tI = S_cI_.T * np.where(S_ > 0, S_tI / S_, 0)[np.newaxis,:]
                R_cS_cI = S_cR_.T * np.where(S_ > 0, (1 - 1/K) * S_cI / S_, 0)[np.newaxis,:]
                R_cS_tI = S_cR_.T * np.where(S_ > 0, S_tI / S_, 0)[np.newaxis,:]

                S_tS_cI = S_tS_ * np.where(S_ > 0, S_cI / S_, 0)[np.newaxis,:]
                S_tS_tI = S_tS_ * np.where(S_ > 0, (1 - 1/(c * N)) * S_tI / S_, 0)[np.newaxis,:]
                I_tS_cI = S_tI_.T * np.where(S_ > 0, S_cI / S_, 0)[np.newaxis,:]
                I_tS_tI = S_tI_.T * np.where(S_ > 0, (1 - 1/(c * N)) * S_tI / S_, 0)[np.newaxis,:]
                R_tS_cI = S_tR_.T * np.where(S_ > 0, S_cI / S_, 0)[np.newaxis,:]
                R_tS_tI = S_tR_.T * np.where(S_ > 0, (1 - 1/(c * N)) * S_tI / S_, 0)[np.newaxis,:]

                
                dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                dc_dt = 2 * np.dot(p, dp_dt)

                dS__dt = sigma * R_ - beta_c * S_cI - beta_t * S_tI
                dI__dt = beta_c * S_cI + beta_t * S_tI - gamma * I_
                dR__dt = gamma * I_ - sigma * R_

                dS_cS__dt = - (beta_c * (S_cS_cI + S_cS_cI.T) + beta_t * (S_cS_tI + S_cS_tI.T)) + sigma * (S_cR_ + S_cR_.T)
                dI_cI__dt = (beta_c * (I_cS_cI + I_cS_cI.T) + beta_t * (I_cS_tI + I_cS_tI.T)) + beta_c * (S_cI_ + S_cI_.T) - 2 * gamma * I_cI_
                dR_cR__dt = gamma * (I_cR_ + I_cR_.T) - 2 * sigma * R_cR_

                dS_cI__dt = (beta_c * (S_cS_cI - I_cS_cI.T) + beta_t * (S_cS_tI - I_cS_tI.T)) - (beta_c + gamma) * S_cI_ + sigma * I_cR_.T
                dS_cR__dt = - (beta_c * R_cS_cI.T + beta_t * R_cS_tI.T) + gamma * S_cI_ + sigma * (R_cR_ - S_cR_)
                dI_cR__dt = (beta_c * R_cS_cI.T + beta_t * R_cS_tI.T) + gamma * I_cI_ - (gamma + sigma) * I_cR_

                dS_tS__dt = - (beta_c * (S_tS_cI + S_tS_cI.T) + beta_t * (S_tS_tI + S_tS_tI.T)) + sigma * (S_tR_ + S_tR_.T) \
                          + dc_dt * np.outer(S_, S_)
                dI_tI__dt = (beta_c * (I_tS_cI + I_tS_cI.T) + beta_t * (I_tS_tI + I_tS_tI.T)) + beta_t * (S_tI_ + S_tI_.T) - 2 * gamma * I_tI_ \
                          + dc_dt * np.outer(I_, I_)
                dR_tR__dt = gamma * (I_tR_ + I_tR_.T) - 2 * sigma * R_tR_ \
                          + dc_dt * np.outer(R_, R_)
                dS_tI__dt = (beta_c * (S_tS_cI - I_tS_cI.T) + beta_t * (S_tS_tI - I_tS_tI.T)) - (beta_t + gamma) * S_tI_ + sigma * I_tR_.T \
                          + dc_dt * np.outer(S_, I_)
                dS_tR__dt = - (beta_c * R_tS_cI.T + beta_t * R_tS_tI.T) + gamma * S_tI_ + sigma * (R_tR_ - S_tR_) \
                          + dc_dt * np.outer(S_, R_)
                dI_tR__dt = (beta_c * R_tS_cI.T + beta_t * R_tS_tI.T) + gamma * I_tI_ - (gamma + sigma) * I_tR_ \
                          + dc_dt * np.outer(I_, R_)

                __dV_dt = np.concatenate((dp_dt, np.vstack([dS__dt, dI__dt, dR__dt, dS_cS__dt, dI_cI__dt, dR_cR__dt, dS_cI__dt, dS_cR__dt, dI_cR__dt, dS_tS__dt, dI_tI__dt, dR_tR__dt, dS_tI__dt, dS_tR__dt, dI_tR__dt]).reshape(-1)))

                return __dV_dt

            V0 = np.concatenate((p0, np.vstack([S_0, I_0, R_0, S_cS_0, I_cI_0, R_cR_0, S_cI_0, S_cR_0, I_cR_0, S_tS_0, I_tI_0, R_tR_0, S_tI_0, S_tR_0, I_tR_0]).reshape(-1)))
            solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degrees), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

            t = solve_ivp_solution.t
            V_t = solve_ivp_solution.y
            
            c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
            M = len(self.community_network_degrees)
            S__t, I__t, R__t, S_cS__t, I_cI__t, R_cR__t, S_cI__t, S_cR__t, I_cR__t, S_tS__t, I_tI__t, R_tR__t, S_tI__t, S_tR__t, I_tR__t, _ = np.vsplit(V_t[self.world_size:,:].reshape(1 * 3 + M * 2*6,M,len(T)), np.cumsum([1] * 3 + [M] * 2*6))
            S__t = S__t.squeeze()
            I__t = I__t.squeeze()
            R__t = R__t.squeeze()

            if normalised_output:
                S__t /= self.community_size
                I__t /= self.community_size
                R__t /= self.community_size

                S_cS__t /= self.community_size * (self.community_size - 1)
                I_cI__t /= self.community_size * (self.community_size - 1)
                R_cR__t /= self.community_size * (self.community_size - 1)
                S_cI__t /= self.community_size * (self.community_size - 1)
                S_cR__t /= self.community_size * (self.community_size - 1)
                I_cR__t /= self.community_size * (self.community_size - 1)
                S_tS__t /= self.community_size * (self.community_size - 1)
                I_tI__t /= self.community_size * (self.community_size - 1)
                R_tR__t /= self.community_size * (self.community_size - 1)
                S_tI__t /= self.community_size * (self.community_size - 1)
                S_tR__t /= self.community_size * (self.community_size - 1)
                I_tR__t /= self.community_size * (self.community_size - 1)

            if full_output:
                return pd.DataFrame(data=dict(zip(('S', 'I', 'R', 'ScS', 'IcI', 'RcR', 'ScI', 'ScR', 'IcR', 'StS', 'ItI', 'RtR', 'StI', 'StR', 'ItR',), (S__t.sum(axis=0), I__t.sum(axis=0), R__t.sum(axis=0), S_cS__t.sum(axis=(0,1)), I_cI__t.sum(axis=(0,1)), R_cR__t.sum(axis=(0,1)), S_cI__t.sum(axis=(0,1)), S_cR__t.sum(axis=(0,1)), I_cR__t.sum(axis=(0,1)), S_tS__t.sum(axis=(0,1)), I_tI__t.sum(axis=(0,1)), R_tR__t.sum(axis=(0,1)), S_tI__t.sum(axis=(0,1)), S_tR__t.sum(axis=(0,1)), I_tR__t.sum(axis=(0,1)),))), index=t)
            else:
                return pd.DataFrame(data=dict(zip(('I',), (I__t.sum(axis=0),))), index=t)


    def SIS_irregular_macroscopic_mean_field_timeseries(self, T, order=1, normalised_output=True, full_output=True):

        # INITIAL CONDITIONS

        if (0 <= self.initial_site < self.world_size):
            p0 = np.zeros(self.world_size)
            p0[self.initial_site] = 1
        else:
            p0 = np.ones(self.world_size) / self.world_size


        c0 = np.dot(p0, p0)
        S_0 = (1 - self.initial_prevalence) * self.community_network_degree_counts
        I_0 = self.initial_prevalence * self.community_network_degree_counts

        S_cI_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(S_0, I_0)
        S_cS_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(S_0, S_0)
        I_cI_0 = np.outer(self.community_network_degrees, self.community_network_degrees) / np.dot(self.community_network_degrees, self.community_network_degree_counts) * np.outer(I_0, I_0)
        
        S_tI_0 = c0 * np.outer(S_0, I_0)
        S_tS_0 = c0 * np.outer(S_0, S_0)
        I_tI_0 = c0 * np.outer(I_0, I_0)



        if order == 1:

            def dV_dt(t, V, mu, beta_c, beta_t, gamma, K):

                p = V[0:self.world_size]

                M = len(K)
                S_, I_, _ = np.vsplit(V[self.world_size:].reshape(1 * 2,M), np.cumsum([1] * 2))
    
                c = np.dot(p, p)
    
                S_cI = (np.outer(K, K) / np.dot(K, self.community_network_degree_counts) * np.outer(S_, I_)).sum(axis=1)
                S_tI = (c * np.outer(S_, I_)).sum(axis=1)


                dp_dt = - mu * self.transport_network_laplacian(t).T @ p

                dS__dt = -(beta_c * S_cI + beta_t * S_tI) + gamma * I_
                dI__dt = (beta_c * S_cI + beta_t * S_tI) - gamma * I_
    
                __dV_dt = np.concatenate((dp_dt, np.vstack([dS__dt, dI__dt]).reshape(-1)))
    
                return __dV_dt

            V0 = np.concatenate((p0, np.vstack([S_0, I_0]).reshape(-1)))
            solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_network_degrees), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)
    
            t = solve_ivp_solution.t
            V_t = solve_ivp_solution.y

            c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
            M = len(self.community_network_degrees)
            S__t, I__t, _ = np.vsplit(V_t[self.world_size:,:].reshape(1 * 2,M,len(T)), np.cumsum([1] * 2))
            S__t = S__t.squeeze()
            I__t = I__t.squeeze()

            
            if normalised_output:
                S__t /= self.community_size
                I__t /= self.community_size

            if full_output:
                return pd.DataFrame(data=dict(zip(('S', 'I',), (S__t.sum(axis=0), I__t.sum(axis=0),))), index=t)
            else:
                return pd.DataFrame(data=dict(zip(('I',), (I__t.sum(axis=0),))), index=t)

        elif order == 2:
            
            def dV_dt(t, V, mu, beta_c, beta_t, gamma, N, K):
                
                p = V[0:self.world_size]
                
                M = len(K)
                S_, I_, S_cI_, S_cS_, I_cI_, S_tI_, S_tS_, I_tI_, _ = np.vsplit(V[self.world_size:].reshape(1 * 2 + M * 2*3,M), np.cumsum([1] * 2 + [M] * 2*3))
                S_ = S_.flatten()
                I_ = I_.flatten()


                c = np.dot(p, p)


                S_cI = S_cI_.sum(axis=1)
                S_tI = S_tI_.sum(axis=1)

                S_cS_cI = S_cS_ * np.where(S_ > 0, (1 - 1/K) * S_cI / S_, 0)[np.newaxis,:]
                S_cS_tI = S_cS_ * np.where(S_ > 0, S_tI / S_, 0)[np.newaxis,:]
                I_cS_cI = S_cI_.T * np.where(S_ > 0, (1 - 1/K) * S_cI / S_, 0)[np.newaxis,:]
                I_cS_tI = S_cI_.T * np.where(S_ > 0, S_tI / S_, 0)[np.newaxis,:]

                S_tS_cI = S_tS_ * np.where(S_ > 0, S_cI / S_, 0)[np.newaxis,:]
                S_tS_tI = S_tS_ * np.where(S_ > 0, (1 - 1/(c * N)) * S_tI / S_, 0)[np.newaxis,:]
                I_tS_cI = S_tI_.T * np.where(S_ > 0, S_cI / S_, 0)[np.newaxis,:]
                I_tS_tI = S_tI_.T * np.where(S_ > 0, (1 - 1/(c * N)) * S_tI / S_, 0)[np.newaxis,:]

                
                dp_dt = - mu * self.transport_network_laplacian(t).T @ p
                dc_dt = 2 * np.dot(p, dp_dt)

                dS__dt = gamma * I_ - beta_c * S_cI - beta_t * S_tI
                dI__dt = beta_c * S_cI + beta_t * S_tI - gamma * I_

                dS_cS__dt = - (beta_c * (S_cS_cI + S_cS_cI.T) + beta_t * (S_cS_tI + S_cS_tI.T)) + gamma * (S_cI_ + S_cI_.T)
                dI_cI__dt = (beta_c * (I_cS_cI + I_cS_cI.T) + beta_t * (I_cS_tI + I_cS_tI.T)) + beta_c * (S_cI_ + S_cI_.T) - 2 * gamma * I_cI_

                dS_cI__dt = (beta_c * (S_cS_cI - I_cS_cI.T) + beta_t * (S_cS_tI - I_cS_tI.T)) - (beta_c + gamma) * S_cI_ + gamma * I_cI_

                dS_tS__dt = - (beta_c * (S_tS_cI + S_tS_cI.T) + beta_t * (S_tS_tI + S_tS_tI.T)) + gamma * (S_tI_ + S_tI_.T) \
                          + dc_dt * np.outer(S_, S_)
                dI_tI__dt = (beta_c * (I_tS_cI + I_tS_cI.T) + beta_t * (I_tS_tI + I_tS_tI.T)) + beta_t * (S_tI_ + S_tI_.T) - 2 * gamma * I_tI_ \
                          + dc_dt * np.outer(I_, I_)
                dS_tI__dt = (beta_c * (S_tS_cI - I_tS_cI.T) + beta_t * (S_tS_tI - I_tS_tI.T)) - (beta_t + gamma) * S_tI_ + gamma * I_tI_.T \
                          + dc_dt * np.outer(S_, I_)

                __dV_dt = np.concatenate((dp_dt, np.vstack([dS__dt, dI__dt, dS_cI__dt, dS_cS__dt, dI_cI__dt, dS_tI__dt, dS_tS__dt, dI_tI__dt]).reshape(-1)))

                return __dV_dt

            V0 = np.concatenate((p0, np.vstack([S_0, I_0, S_cI_0, S_cS_0, I_cI_0, S_tI_0, S_tS_0, I_tI_0]).reshape(-1)))
            solve_ivp_solution = scipy.integrate.solve_ivp(dV_dt, (T.min(), T.max()), V0, args=(self.mobility_rate, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_size, self.community_network_degrees), t_eval=T, method='LSODA', atol=1e-9, rtol=1e-6)

            t = solve_ivp_solution.t
            V_t = solve_ivp_solution.y
            
            c_t = np.sum(V_t[:self.world_size,:]**2, axis=0)
            M = len(self.community_network_degrees)
            S__t, I__t, S_cI__t, S_cS__t, I_cI__t, S_tI__t, S_tS__t, I_tI__t, _ = np.vsplit(V_t[self.world_size:,:].reshape(1 * 2 + M * 2*3,M,len(T)), np.cumsum([1] * 2 + [M] * 2*3))
            S__t = S__t.squeeze()
            I__t = I__t.squeeze()

            if normalised_output:
                S__t /= self.community_size
                I__t /= self.community_size

                S_cI__t /= self.community_size * (self.community_size - 1)
                S_cS__t /= self.community_size * (self.community_size - 1)
                I_cI__t /= self.community_size * (self.community_size - 1)
                S_tI__t /= self.community_size * (self.community_size - 1)
                S_tS__t /= self.community_size * (self.community_size - 1)
                I_tI__t /= self.community_size * (self.community_size - 1)

            if full_output:
                return pd.DataFrame(data=dict(zip(('S', 'I', 'ScI', 'ScS', 'IcI', 'StI', 'StS', 'ItI',), (S__t.sum(axis=0), I__t.sum(axis=0), S_cI__t.sum(axis=(0,1)), S_cS__t.sum(axis=(0,1)), I_cI__t.sum(axis=(0,1)), S_tI__t.sum(axis=(0,1)), S_tS__t.sum(axis=(0,1)), I_tI__t.sum(axis=(0,1)),))), index=t)
            else:
                return pd.DataFrame(data=dict(zip(('I',), (I__t.sum(axis=0),))), index=t)


    def macroscopic_mean_field_endemic_equilibrium(self, order=1, normalised_output=True):
        #print('>> macroscopic_mean_field_endemic_equilibrium')

        self.epidemic_flavour = self.determine_epidemic_flavour()



        if len(self.community_network_degrees) > 1:

            return None

        else:

            ## THIS IS A MEAN FIELD MODEL VALID ONLY FOR k-REGULAR (EPIDEMIC) COMMUNITY NETWORKS
            #
            #k = np.mean([v for _, v in self.community_network.degree()])
            #if np.std([v for _, v in self.community_network.degree()]) > 0:
            #    warnings.warn('The implemented macroscopic mean field model is not valid for non-regular community networks.', UserWarning)
            #
            #assert all([_ == k for _ in self.community_network_degree])


            if self.epidemic_flavour in ['SIRS', 'SI', 'SIR']:
                return self.SIRS_macroscopic_mean_field_endemic_equilibrium(order=order, normalised_output=normalised_output)

            if self.epidemic_flavour == 'SIS':
                return self.SIS_macroscopic_mean_field_endemic_equilibrium(order=order, normalised_output=normalised_output)


        raise NotImplementedError('An {} mean-field model has not been implemented yet.'.format(self.epidemic_flavour))


    def SIRS_macroscopic_mean_field_endemic_equilibrium(self, order=1, normalised_output=True):

        # TRANSPORT EQUILIBRIUM DISTRIBUTION

        p_inf = self.transport_network_equilibrium_distributions[0]



        if order == 1:
            
            c_inf = np.dot(p_inf, p_inf)
            X_inf = self.community_size * (self.community_infection_rate * self.community_network_degree[0] / self.community_size + self.transport_infection_rate * c_inf) / self.recovery_rate

            if X_inf > 1:
                S_inf = self.community_size * 1 / X_inf
                I_inf = self.community_size * self.immunity_loss_rate / (self.recovery_rate + self.immunity_loss_rate) * (1 - 1 / X_inf)
                R_inf = self.community_size * self.recovery_rate / (self.recovery_rate + self.immunity_loss_rate) * (1 - 1 / X_inf)
            else:
                S_inf = self.community_size
                I_inf = 0
                R_inf = 0

            if normalised_output:
                S_inf /= self.community_size
                I_inf /= self.community_size
                R_inf /= self.community_size

            return {'p': p_inf, 'S' : S_inf, 'I' : I_inf, 'R' : R_inf}

        elif order == 2:

            # INITIAL CONDITIONS

            c0 = np.dot(p_inf, p_inf)

            _ = self.SIRS_macroscopic_mean_field_endemic_equilibrium(order=1, normalised_output=False)

            S0 = _['S']
            I0 = _['I']
            R0 = _['R']

            ScS0 = self.community_network_degree[0] / self.community_size * S0**2
            IcI0 = self.community_network_degree[0] / self.community_size * I0**2
            RcR0 = self.community_network_degree[0] / self.community_size * R0**2
            ScI0 = self.community_network_degree[0] / self.community_size * S0 * I0
            ScR0 = self.community_network_degree[0] / self.community_size * S0 * R0
            IcR0 = self.community_network_degree[0] / self.community_size * I0 * R0

            StS0 = c0 * S0**2
            ItI0 = c0 * I0**2
            RtR0 = c0 * R0**2
            StI0 = c0 * S0 * I0
            StR0 = c0 * S0 * R0
            ItR0 = c0 * I0 * R0

            
            def dF_dt(X, p_inf, beta_c, beta_t, gamma, sigma, N, k):
    
                if gamma == 0:
                    sigma = 0

                S, I, ScS, IcI, ScI, StS, ItI, StI = X

                c_inf = np.dot(p_inf, p_inf)

                R = N - (S + I)
                ScR = k * S - ScS - ScI
                IcR = k * I - IcI - ScI
                RcR = k * R - ScR - IcR
                StR = c_inf * N * S - StS - StI
                ItR = c_inf * N * I - ItI - StI
                RtR = c_inf * N * R - StR - ItR


                if S > 0:
                    ScScI = (1 - 1/k) * ScS * ScI / S
                    ScStI = ScS * StI / S
                    IcScI = (1 - 1/k) * ScI**2 / S
                    IcStI = ScI * StI / S
                    RcScI = (1 - 1/k) * ScR * ScI / S
                    RcStI = ScR * StI / S

                    StScI = StS * ScI / S
                    StStI = (1 - 1/(c_inf * N)) * StS * StI / S
                    ItScI = StI * ScI / S
                    ItStI = (1 - 1/(c_inf * N)) * StI**2 / S
                    RtScI = StR * ScI / S
                    RtStI = (1 - 1/(c_inf * N)) * StR * StI / S
                else:
                    ScScI = 0
                    ScStI = 0
                    IcScI = 0
                    IcStI = 0
                    RcScI = 0
                    RcStI = 0

                    StScI = 0
                    StStI = 0
                    ItScI = 0
                    ItStI = 0
                    RtScI = 0
                    RtStI = 0

                    
                dS_dt = sigma * R - beta_c * ScI - beta_t * StI
                dI_dt = beta_c * ScI + beta_t * StI - gamma * I

                dScS_dt = -2 * (beta_c * ScScI + beta_t * ScStI) + 2 * sigma * ScR
                dIcI_dt = 2 * (beta_c * IcScI + beta_t * IcStI) + 2 * beta_c * ScI - 2 * gamma * IcI
                dScI_dt = (beta_c * (ScScI - IcScI) + beta_t * (ScStI - IcStI)) - (beta_c + gamma) * ScI + sigma * IcR

                dStS_dt = -2 * (beta_c * StScI + beta_t * StStI) + 2 * sigma * StR
                dItI_dt = 2 * (beta_c * ItScI + beta_t * ItStI) + 2 * beta_t * StI - 2 * gamma * ItI
                dStI_dt = (beta_c * (StScI - ItScI) + beta_t * (StStI - ItStI)) - (beta_t + gamma) * StI + sigma * ItR

                return np.array([dS_dt, dI_dt, dScS_dt, dIcI_dt, dScI_dt, dStS_dt, dItI_dt, dStI_dt])


            root_finder = scipy.optimize.root(dF_dt, np.array([S0, I0, ScS0, IcI0, ScI0, StS0, ItI0, StI0]), args=(p_inf, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.immunity_loss_rate, self.community_size, self.community_network_degree[0]))


            if root_finder.success:
                S_inf, I_inf, ScS_inf, IcI_inf, ScI_inf, StS_inf, ItI_inf, StI_inf = root_finder.x
            else:
                warnings.warn('2nd-order mean-field equilibrium solution could not be found.')
                S_inf, I_inf, ScS_inf, IcI_inf, ScI_inf, StS_inf, ItI_inf, StI_inf = np.array([S0, I0, ScS0, IcI0, ScI0, StS0, ItI0, StI0])
                
            c_inf = np.dot(p_inf, p_inf)

            R_inf = self.community_size - (S_inf + I_inf)
            ScR_inf = self.community_network_degree[0] * S_inf - ScS_inf - ScI_inf
            IcR_inf = self.community_network_degree[0] * I_inf - IcI_inf - ScI_inf
            RcR_inf = self.community_network_degree[0] * R_inf - ScR_inf - IcR_inf
            StR_inf = c_inf * self.community_size * S_inf - StS_inf - StI_inf
            ItR_inf = c_inf * self.community_size * I_inf - ItI_inf - StI_inf
            RtR_inf = c_inf * self.community_size * R_inf - StR_inf - ItR_inf

            if normalised_output:
                S_inf /= self.community_size
                I_inf /= self.community_size
                R_inf /= self.community_size

                ScS_inf /= self.community_size * (self.community_size - 1)
                IcI_inf /= self.community_size * (self.community_size - 1)
                RcR_inf /= self.community_size * (self.community_size - 1)
                ScI_inf /= self.community_size * (self.community_size - 1)
                ScR_inf /= self.community_size * (self.community_size - 1)
                IcR_inf /= self.community_size * (self.community_size - 1)
                StS_inf /= self.community_size * (self.community_size - 1)
                ItI_inf /= self.community_size * (self.community_size - 1)
                RtR_inf /= self.community_size * (self.community_size - 1)
                StI_inf /= self.community_size * (self.community_size - 1)
                StR_inf /= self.community_size * (self.community_size - 1)
                ItR_inf /= self.community_size * (self.community_size - 1)

            return {'p': p_inf, 'S' : S_inf, 'I' : I_inf, 'R' : R_inf, 'ScS' : ScS_inf, 'IcI' : IcI_inf, 'RcR' : RcR_inf, 'ScI' : ScI_inf, 'ScR' : ScR_inf, 'IcR' : IcR_inf, 'StS' : StS_inf, 'ItI' : ItI_inf, 'RtR' : RtR_inf, 'StI' : StI_inf, 'StR' : StR_inf, 'ItR' : ItR_inf}
            

    def SIS_macroscopic_mean_field_endemic_equilibrium(self, order=1, normalised_output=True):

        # TRANSPORT EQUILIBRIUM DISTRIBUTION

        p_inf = self.transport_network_equilibrium_distributions[0]



        if order == 1:
            
            c_inf = np.dot(p_inf, p_inf)
            X_inf = self.community_size * (self.community_infection_rate * self.community_network_degree[0] / self.community_size + self.transport_infection_rate * c_inf) / self.recovery_rate

            if X_inf > 1:
                S_inf = self.community_size * 1 / X_inf
                I_inf = self.community_size * (1 - 1 / X_inf)
            else:
                S_inf = self.community_size
                I_inf = 0

            if normalised_output:
                S_inf /= self.community_size
                I_inf /= self.community_size

            return {'p': p_inf, 'S' : S_inf, 'I' : I_inf}

        elif order == 2:

            # INITIAL CONDITIONS

            c0 = np.dot(p_inf, p_inf)

            _ = self.SIS_macroscopic_mean_field_endemic_equilibrium(order=1, normalised_output=False)

            S0 = _['S']
            I0 = _['I']

            ScS0 = self.community_network_degree[0] / self.community_size * S0**2
            IcI0 = self.community_network_degree[0] / self.community_size * I0**2
            ScI0 = self.community_network_degree[0] / self.community_size * S0 * I0

            StS0 = c0 * S0**2
            ItI0 = c0 * I0**2
            StI0 = c0 * S0 * I0

            
            def dF_dt(X, p_inf, beta_c, beta_t, gamma, N, k):

                S, ScS, StS = X

                c_inf = np.dot(p_inf, p_inf)

                I = N - S
                ScI = k * S - ScS
                IcI = k * I - ScI
                StI = c_inf * N * S - StS
                ItI = c_inf * N * I - StI


                if S > 0:
                    ScScI = (1 - 1/k) * ScS * ScI / S
                    StScI = StS * ScI / S
                    ScStI = ScS * StI / S
                    StStI = (1 - 1/(c_inf * N)) * StS * StI / S
                    IcScI = (1 - 1/k) * ScI**2 / S
                    IcStI = ScI * StI / S
                    ItStI = (1 - 1/(c_inf * N)) * StI**2 / S
                else:
                    ScScI = 0
                    StScI = 0
                    ScStI = 0
                    StStI = 0
                    IcScI = 0
                    IcStI = 0
                    ItStI = 0
                    

                dS_dt = gamma * I - beta_c * ScI - beta_t * StI

                dScS_dt = -2 * (beta_c * ScScI + beta_t * ScStI) + 2 * gamma * ScI

                dStS_dt = -2 * (beta_c * StScI + beta_t * StStI) + 2 * gamma * StI

                return np.array([dS_dt, dScS_dt, dStS_dt])


            root_finder = scipy.optimize.root(dF_dt, np.array([S0, ScS0, StS0]), args=(p_inf, self.community_infection_rate, self.transport_infection_rate, self.recovery_rate, self.community_size, self.community_network_degree[0]))


            if root_finder.success:
                S_inf, ScS_inf, StS_inf = root_finder.x
            else:
                warnings.warn('2nd-order mean-field equilibrium solution could not be found.')
                S_inf, ScS_inf, StS_inf = np.array([S0, ScS0, StS0])
                
            c_inf = np.dot(p_inf, p_inf)

            I_inf = self.community_size - S_inf
            ScI_inf = self.community_network_degree[0] * S_inf - ScS_inf
            IcI_inf = self.community_network_degree[0] * I_inf - ScI_inf
            StI_inf = c_inf * self.community_size * S_inf - StS_inf
            ItI_inf = c_inf * self.community_size * I_inf - StI_inf

            if normalised_output:
                S_inf /= self.community_size
                I_inf /= self.community_size

                ScS_inf /= self.community_size * (self.community_size - 1)
                IcI_inf /= self.community_size * (self.community_size - 1)
                ScI_inf /= self.community_size * (self.community_size - 1)
                StS_inf /= self.community_size * (self.community_size - 1)
                ItI_inf /= self.community_size * (self.community_size - 1)
                StI_inf /= self.community_size * (self.community_size - 1)

            return {'p': p_inf, 'S' : S_inf, 'I' : I_inf, 'ScS' : ScS_inf, 'IcI' : IcI_inf, 'ScI' : ScI_inf, 'StS' : StS_inf, 'ItI' : ItI_inf, 'StI' : StI_inf}
