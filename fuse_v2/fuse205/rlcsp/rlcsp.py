#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 10:41:17 2020
CSP variation of Reinforce algorithm designed to work with basin-hopping CSP codes
on the step where the next action to change the trial structure is being selected
@author: Elena Zamaraeva
"""
import numpy as np
import math
import statistics
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import os
from .state import State
from matplotlib.animation import FuncAnimation, PillowWriter
import time
#from mysql.connector import connect, Error
import sqlite3
import random
import string
import scipy.stats
import heapq

"""
    This file contains a Reinforce class used to learn the dynamic policy
"""


class Reinforce:

    # Three possible types of rewards
    unscaled = 'unscaled'  # unscaled energy drop
    change_in_features = 'features'  # energy drop as a change in features
    uniqueness = 'uniqueness'  # binary reward for uniqueness of the energy

    f_energy = 'energy'
    f_sl = 'sl'
    f_old_energy = 'old_energy'

    def __init__(self, actions=["T1", "T2", "T3", "T4", "T5", "T6"],
                 alpha=0,  # the learning rate
                 params_db={'host': 'localhost',
                            'database': '',
                            'user': '',
                            'password': ''},
                 reinforce_table='reinforce',  # the table containing features history
                 theta_table='theta',  # the table containing theta history
                 reward_type='',  # the reward type
                 sl_range=[5, 10],  # possible stack lengths (used for MC-EMMA)
                 features_set=['energy'],  # available features
                 episode_length=1,  # episode length used if episodes have fixed length
                 fixed_episodes=True,  # indicates if episodes have fixed length
                 max_energy=0,  # maximum allowed energy
                 reinforce_id=0,  # agent ID used for the policy sharing
                 db_max_attempts=5,  # maximum number of attempts to send SQL queries
                 debug=False,  # indicates if the agent works in degubbing mode and print more information
                 reg_params={'e_threshold': 0,  # the entropy threshold below which entropy regularization turns on
                             'beta': 0,  # the beta parameter in the entropy regularization mechanizm
                             'free_term': 0,  # free term in h function
                             'h_type': 'non-linear',  # the type of h function
                             'theta': '',  # policy vector can be used to define the starting policy to be non-uniform
                             'zero_reward_penalty': 0,  # the penalty for the zero reward
                             'non_converge_penalty': 0,  # the penalty for the non-converged energy
                             'non_unique_penalty': 0,  # the penalty for the non-unique structure
                             'scale_reward': False,  # indicates whether the reward is being scaled
                             'epsilon': 0,  # epsilon parameter for the epsilon-greedy policy
                             'non_unique_reward': True,  # indicates if the reward for non_unique energy is summarized
                                                         # with the non-unique penalty
                             'step_reward_limit': 100,  # the maximum number of rewards used for the scaling
                             'reward_limit_last': False,  # if True we take the _last_ step_reward_limit rewards
                                                            # for the reward scaling
                             'smart_penalty': False, # indicates if penalties are calculated by the agent
                             'storage_type': 'mysql'}):  # "mysql" or "sqlite" for different storage types

        self.reward_type = Reinforce.change_in_features if reward_type == '' else reward_type
        if 'zero_reward_penalty' not in reg_params:
            self.zero_reward_penalty = 0
        else:
            self.zero_reward_penalty = reg_params['zero_reward_penalty']
        if 'non_converge_penalty' not in reg_params:
            self.non_converge_penalty = 0
        else:
            self.non_converge_penalty = reg_params['non_converge_penalty']
        if 'non_unique_penalty' not in reg_params:
            self.non_unique_penalty = 0
        else:
            self.non_unique_penalty = reg_params['non_unique_penalty']

        self.max_energy = max_energy

        self.beta = 1

        # regularization parameters: e_threshold is the entropy threshold after which we apply regularization
        # beta is the coefficient of confidence penalty
        self.reg_params = reg_params
        if 'e_threshold' not in self.reg_params:
            self.reg_params['e_threshold'] = 0
        if 'epsilon' not in self.reg_params:
            self.reg_params['epsilon'] = 0
        if 'beta' not in self.reg_params:
            self.reg_params['beta'] = 0
        if 'scale_reward' not in self.reg_params:
            self.reg_params['scale_reward'] = False
        if 'non_unique_reward' not in self.reg_params:
            self.reg_params['non_unique_reward'] = True
        if 'reward_limit_last' not in self.reg_params:
            self.reg_params['reward_limit_last'] = False
        if 'smart_penalty' not in self.reg_params:
            self.reg_params['smart_penalty'] = False
        
        self.reg_params['storage_type'] = self.reg_params.get("storage_type", 'mysql')

        self.free_term = self.reg_params['free_term']
        self.h_type = self.reg_params['h_type']

        self.episode_length = episode_length
        self.fixed_episodes = fixed_episodes
        self.episode = []

        self.id = reinforce_id
        self.db_max_attempts = db_max_attempts

        self.actions = actions
        self.theta = [0 for k in range(len(features_set) * len(self.actions) * 2)]

        self.thetas = []
        self.f_array = {f: [] for f in features_set}
        self.f_mean = {f: None for f in features_set}
        self.f_std = {f: None for f in features_set}
        self.f_min = {f: None for f in features_set}
        self.f_max = {f: None for f in features_set}

        self.mean_reward = None
        self.std_reward = None

        self.debug = debug

        if Reinforce.f_sl in features_set:
            self.f_min[Reinforce.f_sl] = min(sl_range)
            self.f_max[Reinforce.f_sl] = max(sl_range)

        self.rewards = []
        self.alpha = alpha
        self.weight = np.zeros(1)

        self.step = 0

        # maximum steps for updating mean_energy. After this number of states we stop updating self.mean_energy
        self.step_energy_limit = 10000

        if 'theta' in self.reg_params:
            if self.reg_params['theta'] != '':
                self.theta = self.reg_params['theta']

        # create and fill DB with theta if it does not exist
        random_db_name = False
        if 'database' in params_db:
            if params_db['database'] == '':
                random_db_name = True
        else:
            random_db_name = True

        if random_db_name:
            letters = string.ascii_lowercase
            params_db['database'] = 'reinforce_' + ''.join(random.choice(letters) for i in range(5))

        self.params_db = params_db
        self.theta_table = theta_table
        self.reinforce_table = reinforce_table
        self.conn_db = None

        self.init_sql_queries()

        self.fill_params_db()

    def init_sql_queries(self):
        """
        Fills in the SQL queries templates
        """
        self.sql_select_theta = 'SELECT theta FROM ' + self.theta_table + ' ORDER BY id DESC LIMIT 1'
        self.sql_select_theta_history = 'SELECT theta FROM ' + self.theta_table

        self.sql_select_energy_history = 'SELECT new_energy FROM ' + self.reinforce_table + ' LIMIT '

        if self.reg_params['reward_limit_last']:
            self.sql_select_reward_history = 'SELECT unscaled_reward FROM ' + self.theta_table + ' ORDER BY id DESC LIMIT '
        else:
            self.sql_select_reward_history = 'SELECT unscaled_reward FROM ' + self.theta_table + ' LIMIT '

    def sql_update_theta(self, unscaled_reward=0, scaled_reward=0, opt='True'):
        """
        Inserts new policy parameter (theta) in DB
        """
        if self.reg_params['scale_reward']:
            return f"""INSERT INTO {self.theta_table}(theta, reinforce_id, unscaled_reward, scaled_reward, opt) 
                    VALUES('{json.dumps(self.theta)}', {self.id}, {unscaled_reward}, {scaled_reward}, {opt})"""

        return f'INSERT INTO {self.theta_table}(theta, reinforce_id) VALUES("{json.dumps(self.theta)}", {self.id})'

    def sql_update_reinforce(self, old_state, new_state, action, alpha):
        """
        Inserts new state (features parameterized) in DB
        """

        f_tablenames, values = self.features_to_tablenames()
        separator = ', '
        features = separator.join(f_tablenames)
        for i in range(len(values)):
            values[i] = '' if values[i] is None else str(values[i])
        values = separator.join(values)
        json_old_struct = old_state.struct_to_str()
        json_new_struct = new_state.struct_to_str()
        return f'''INSERT INTO {self.reinforce_table}(step,
                                                    old_energy, 
                                                    new_energy, 
                                                    action, 
                                                    step_size,
                                                    old_struct,
                                                    new_struct,
                                                    {features}) VALUES 
                                                    ({self.step},
                                                    {old_state.energy},
                                                    {new_state.energy},
                                                    "{action}",
                                                    {alpha},
                                                    "{json_old_struct}",
                                                    "{json_new_struct}", 
                                                    {values})'''

    def features_to_tablenames(self):
        """
        Translates features names to the corresponding tables names
        """

        tablenames = []
        values = []

        features = sorted([*self.f_array])

        for f in features:
            if f in self.f_max.keys():
                tablenames.append('max_' + f)
                values.append(self.f_max[f])
            if f in self.f_min.keys():
                tablenames.append('min_' + f)
                values.append(self.f_min[f])
            if f in self.f_mean.keys():
                tablenames.append('mean_' + f)
                if self.f_mean[f] is not None:
                    values.append(self.f_mean[f])
            if f in self.f_std.keys():
                tablenames.append('std_' + f)
                values.append(self.f_std[f])

        return tablenames, values

    def select_action(self, state, excluded_actions=[]):
        """
        Selects action according to the given state and current policy excluding all actions from
        the excluded_actions list
        """
        self.load_params()
        action_space = np.array(self.actions)
        permitted_actions = [x for x in self.actions if x not in excluded_actions]
        # if the list of excluded actions contains all actions we permit all actions
        if len(permitted_actions) == 0:
            permitted_actions = self.actions
            excluded_actions = []

        # first return randomly chosen action if state is not valid (energy is > max)
        # or with epsilon probability
        if state.energy > self.max_energy or self.reg_params['epsilon'] > random.random():
            action_space = np.array(permitted_actions)
            return np.random.choice(action_space)

        # if energy is <= max, then calculate the probabilities and select an action
        action_probs = [0 for k in self.actions]
        for i in range(len(action_probs)):
            action_probs[i] = self.action_prob(self.actions[i], state, excluded_actions)
        try:
            action = np.random.choice(action_space,
                                      p=action_probs)
        except ValueError as e:
            print(f'Probabilities {action_probs} with error: {e}')
            action = np.random.choice(action_space)

        return action

    def calc_f_scaling(self):
        """
        Updates features scaling parameters: mean, std, min, max
        """

        features = [*self.f_array]
        for f in features:
            if f is Reinforce.f_energy or f is Reinforce.f_old_energy:
                # if len(self.f_array[f]) == 0:
                if len(self.f_array[f]) > 1:
                    self.f_mean[f] = statistics.mean(self.f_array[f])
                    self.f_std[f] = np.std(self.f_array[f])
                    self.f_max[f] = max(self.f_array[f])
                    self.f_min[f] = min(self.f_array[f])

    def update_f_scaling(self, old_state, new_state, action, alpha):
        """
        Reads and scales features from the new state and saves the state to DB
        """

        if self.step < self.step_energy_limit:
            features = [*self.f_array]
            for f in features:
                if f is Reinforce.f_energy or f is Reinforce.f_old_energy:
                    # if len(self.f_array[f]) == 0:
                    if self.f_mean[f] is None:
                        if self.debug:
                             print('update_f_scaling: Mean is None')
                        self.f_array[f] = [old_state.energy]

                    if new_state.energy < self.max_energy:
                        self.f_array[f].append(new_state.energy)

            self.calc_f_scaling()

        if new_state.energy < self.max_energy:
            self.sql_execute(self.sql_update_reinforce(old_state, new_state, action, alpha))

    def update(self, action, old_state, new_state, opt,end_episode=False):
        """
        Observes new state and end episode if the corresponding parameter is True
        and episodes have not fixed length or if it's just time to end the episode.
        Updates the policy if the episode ends.
        """

        # if the energy of the structure is valid add it to the episode
        if old_state.energy <= self.max_energy:
            if new_state.energy <= self.max_energy or self.non_converge_penalty != 0:
                self.episode.append((old_state, new_state, action))

        # if the episode has the fixed length and achieves it, we end the episode
        # if the episode has the flexible length and we force to end it, end it
        if (len(self.episode) >= self.episode_length and self.fixed_episodes) \
                or (end_episode is True and not self.fixed_episodes):
            return self.end_episode(opt=opt)

    def end_episode(self, opt):
        """
        Ends the current episode, calculates the reward, and updates the policy
        """

        reward = 0

        # update the policy for each step of episode
        for step in range(len(self.episode)):

            old_state, new_state, action = self.episode[step]

            reward = self.episodic_reward(step, log=True)
            # print(old_state)
            # print(new_state)
            # print(reward)

            pi = self.action_prob(action=action, state=old_state)
            if pi == 0:
                print("policy is zero. do not update policy")
                print("Action is %s, old energy is %f, new energy is %f" % (action,
                                                                            old_state.energy,
                                                                            new_state.energy))
                return reward

            self.load_params()

            unscaled_reward = reward
            if self.reg_params['scale_reward'] and self.std_reward not in (None, 0) and self.mean_reward is not None:
                reward = (unscaled_reward - self.mean_reward) / self.std_reward
                if self.debug:
                      print('Reward standardization: ' + str(unscaled_reward) + ' --> ' + str(reward))

            if 'max_reward' in self.reg_params:
                reward = min(reward, self.reg_params['max_reward'])
            if 'min_reward' in self.reg_params:
                reward = max(reward, self.reg_params['max_reward'])

            if self.alpha != 0:
                alpha = self.alpha
            else:
                alpha = self.alpha_func()

            self.update_f_scaling(old_state, new_state, action, alpha)
            self.step += 1

            if self.debug:
                print('Reward: ' + str(reward))

            # on the first steps we force reward to be from -1 to +1
            # because features can be standardized using mean and std of energy by too few samples
            if self.step < 10:
                reward = 1 if reward > 1 else reward
                reward = -1 if reward < -1 else reward
                if self.debug:
                         print('Reward forced to -1+1: ' + str(reward))

            if self.debug:
                print('Alpha: ' + str(alpha))
                print('Old theta: ' + str(self.theta))

            diff = ''
            for i in range(len(self.theta)):
                diff_action_prob = self.diff_action_prob(action=action, state=old_state, diff_var=i)
                diff += str(diff_action_prob) + ', '

                if diff_action_prob != 0:
                    entr = self.diff_entropy(action=action, state=old_state, diff_var=i)

                    if self.step > 10:
                        if self.debug:
                            print(f'Diff var is {i}')
                            print(f'Whole update: alpha * (reward * diff_action_prob - beta * entr) / pi')
                            print(f'Whole update: {alpha} * '
                                  f'({reward} * {diff_action_prob} + {self.reg_params["beta"]} * {entr} ) / {pi}')

                        self.theta[i] += alpha * (reward * diff_action_prob + self.reg_params['beta'] * entr) / pi

            if self.debug:
                print('reward: ' + str(reward))
                print('beta: ' + str(self.reg_params['beta']))
                print('Diff: ' + diff)
                print('Probability: ' + str(pi))
                print('New theta: ' + str(self.theta))

            # update theta param in DB
            self.sql_execute(self.sql_update_theta(unscaled_reward=unscaled_reward, scaled_reward=reward,opt=opt))
            self.thetas.append(self.theta.copy())

        # refresh episode
        self.episode = []
        return reward

    def diff_entropy(self, state, action, diff_var):
        """
        Calculates the entropy derivative
        """

        probs = list(map(self.action_prob, self.actions, [state for i in self.actions]))
        entropy = scipy.stats.entropy(probs)
        action_prob = self.action_prob(action=action, state=state)
        diff_action_prob = self.diff_action_prob(action, state, diff_var)

        if self.debug:
            if diff_action_prob != 0:
                print(f'New entropy diff:')
                print(f'Formula: -diff_action_prob * (math.log(action_prob) + 1) ')
                print(f'Formula: -{diff_action_prob} * ({math.log(action_prob)} + 1) ')
                print(f'diff_entropy is {-diff_action_prob * (math.log(action_prob) + 1)}')

        if entropy < self.reg_params['e_threshold']:
            try:
                return -diff_action_prob * (math.log(action_prob) + 1)
            except ValueError as e:
                e_text = f'{e}: state={str(state)} action={action} diff_var={diff_var} action_prob={action_prob} ' \
                         f'diff_action_prob={diff_action_prob}'
                print(e_text)
                raise
        else:
            return 0

    def reward(self, old_state, new_state, log=False):
        """
        Calculates the reward using features of old and nuew state and penalties if needed
        """

        if self.debug:
            print('Old state energy: ' + str(old_state.energy))
            print('New state energy: ' + str(new_state.energy))

        if len(self.f_array[Reinforce.f_energy]) == 0:

            return 0

        # the penalty can be calculated by the agent
        if self.reg_params['smart_penalty']:
            if old_state.energy == new_state.energy or new_state.energy > self.max_energy or not new_state.isunique:
                energies_num = max(1, int(len(self.f_array[Reinforce.f_energy]) / 10))
                high_energy = statistics.mean(heapq.nlargest(energies_num, self.f_array[Reinforce.f_energy]))
                fake_state = new_state.copy()
                fake_state.energy = high_energy
                res = self.feature(Reinforce.f_energy, old_state, log=False) \
                      - self.feature(Reinforce.f_energy, fake_state, log=False)
                if self.debug:
                     print(f'High energy {fake_state.energy} ')
                     print(f'Penalty is {res} ')
                return res

        # or penalty can be fixed by user
        if old_state.energy == new_state.energy:
            res = self.zero_reward_penalty

        elif new_state.energy > self.max_energy:
            res = self.non_converge_penalty

        else:
            # reward calculation by features depends on the reward type and can be the change in energy
            if self.reward_type is Reinforce.unscaled:
                res = (old_state.energy - new_state.energy)

            # the change in features
            elif self.reward_type is Reinforce.change_in_features:
                res = self.feature(Reinforce.f_energy, old_state, log=log) - self.feature(Reinforce.f_energy, new_state,
                                                                                          log=log)

            # or just +/- penalty for non-uniqueness
            elif self.reward_type is Reinforce.uniqueness:
                res = 0
                if new_state.isunique:
                    res -= self.non_unique_penalty

            if self.non_unique_penalty != 0 and not new_state.isunique:
                if self.reg_params['non_unique_reward']:
                    res += self.non_unique_penalty
                else:
                    res = self.non_unique_penalty
                if self.debug:
                    print(f'Penalty for non-uniqueness, reward: {res}')

        return res

    def episodic_reward(self, step, log=False):
        """
        Calculates the episodic reward
        """

        res = 0
        discount = 1

        for i in range(step, len(self.episode)):
            old_state, new_state, action = self.episode[i]
            res += (discount ** i) * self.reward(old_state, new_state, log=log)

        return res

    def standardize_feature(self, f, val, default_val=0, log=False):
        """
        Z-score normalization for feature scaling
        """

        if log and self.debug:
            print('Feature standardisation')
        if self.f_mean[f] is not None and self.f_std[f] != 0 and self.f_std[f] is not None:
            if log and self.debug:
                print('Mean energy: ' + str(self.f_mean[f]))
                print('Std energy: ' + str(self.f_std[f]))
                print('Current energy: ' + str(val))
                print('Feature value: ' + str((val - self.f_mean[f]) / self.f_std[f]))
            return (val - self.f_mean[f]) / self.f_std[f]
        else:
            if log and self.debug:
                print('Default energy: ' + str(default_val))
            return default_val

    def normalize_feature(self, f, val, default_val=0, min_val=0, max_val=1):
        """
        Normalization for feature scaling to be from min_val to max_val
        """

        if self.f_max[f] is not None and self.f_min[f] is not None:
            return min_val + (val - self.f_min[f]) * (max_val - min_val) / (self.f_max[f] - self.f_min[f])
        else:
            return default_val

    def feature(self, f, state, log=False):
        """
        Returns scaled feature given the state
        """

        if f is Reinforce.f_energy:
            return self.standardize_feature(f, state.energy, default_val=state.energy, log=log)

        if f is Reinforce.f_old_energy:
            return self.standardize_feature(f, state.old_energy, default_val=state.old_energy, log=log)

        if f is Reinforce.f_sl:
            return self.normalize_feature(f, state.sl(), default_val=state.sl(), min_val=-1, max_val=-1)

        return None

    def features(self, state, action):
        """
        Returns features vector for given state and action ready for the multiplication by the theta vector.
        The vector size is number of actions * number of features
        The non-zero coordinates relate to the given action's index
        """

        features = sorted(self.f_array.keys())

        res = np.zeros(len(features) * len(self.actions))

        f_counter = 0

        # fill energy features
        if Reinforce.f_energy in features:
            f = Reinforce.f_energy
            res[len(features) * self.actions.index(action) + f_counter] = self.standardize_feature(f, state.energy)
            f_counter += 1

        # fill old energy features
        if Reinforce.f_old_energy in features:
            f = Reinforce.f_old_energy
            res[len(features) * self.actions.index(action) + f_counter] = self.standardize_feature(f, state.old_energy)
            f_counter += 1

        # fill stack length features
        if Reinforce.f_sl in features:
            f = Reinforce.f_sl
            res[len(features) * self.actions.index(action) + f_counter] = self.normalize_feature(f, state.sl(),
                                                                                                 default_val=0,
                                                                                                 min_val=-1,
                                                                                                 max_val=+1)
            f_counter += 1
        return res

    def action_prob(self, action, state, excluded_actions=[]):
        """
        Returns probability to be selected given the action and the state, softmax function
        """
        if action in excluded_actions:
            return 0
        try:
            h_list = []
            h = 0
            for a in self.actions:
                if a not in excluded_actions:
                    h_list.append(self.beta * self.h(state, a))
                    if a == action:
                        h = h_list[-1]

            h_list = np.array(h_list)
            max_h = h_list.max()
            h -= max_h
            h_list -= max_h
            return math.exp(h) / np.exp(h_list).sum()

        except OverflowError as e:
            raise Exception(
                f'''action_prob error: {e}, step is {self.step}, h is {self.h(state, action)}, action is {action}, \nstate is {state},
                    Features are {self.features(state, action)}, \ntheta is {self.theta}''')

        return 0

    @staticmethod
    def st_action_prob(action, state, theta, actions,
                       f_mean, f_std, f_max, f_min, f_array,
                       excluded_actions=[]):
        """
        Returns probability to be selected given the action, the state and all other parameters.
        Can be used without the agent
        """

        reinforce = Reinforce(actions=actions, theta=theta, features_set=f_array.keys())
        reinforce.f_array = f_array
        reinforce.f_mean = f_mean
        reinforce.f_std = f_std
        reinforce.f_min = f_min
        reinforce.f_max = f_max
        return reinforce.action_prob(action, state, excluded_actions)

    def exp_sum(self, state, excluded_actions=[]):
        """
        Sum of exponents for softmax function
        """

        res = 0
        for i in range(len(self.actions)):
            if not self.actions[i] in excluded_actions:
                res += math.exp(self.beta * self.h(state, self.actions[i]))

        return res

    def diff_action_prob(self, action, state, diff_var):
        """
        Derivative of policy on diff_var param
        """

        diff_h = -1

        try:
            diff_h = self.diff_h(state, action, diff_var)

            action_prob = self.action_prob(action, state)

            return action_prob * diff_h * (1 - action_prob)

        except OverflowError as e:
            print(f'diff_action_prob function error: {e}, diff_h is {diff_h}')
            raise Exception(f'''diff_action_prob function error: {e}, 
                                diff_h is {diff_h},
                                step is {self.step}, h is {self.h(state, action)}, 
                                action is {action}, \nstate is {state},
                                Features are {self.features(state, action)}, \ntheta is {self.theta}''')
        except ZeroDivisionError as e:
            raise Exception(f'''diff_action_prob function error: {e}, 
                                diff_h is {diff_h},
                                step is {self.step}, h is {self.h(state, action)}, 
                                action is {action}, \nstate is {state},
                                Features are {self.features(state, action)}, \ntheta is {self.theta}''')

    def diff_action_prob_vector(self, action, state):
        """
        Derivative of policy as a vector
        """

        diff_vector = []

        for diff_var in range(len(self.theta)):
            diff_vector.append(self.diff_action_prob(action, state, diff_var))

        return diff_vector

    def h(self, state, action):
        """
        Action preferences function
        """
        features = self.features(state, action)
        start_coef = []
        for i in range(len(features)):
            if features[i] != 0:
                start_coef.append(self.theta[2 * i + 1])
                if self.h_type != 'linear':
                    features[i] = np.sign(features[i]) * (np.abs(features[i])) ** (1 / 3)
        return np.multiply(self.theta[::2], features).sum() + sum(start_coef) * self.free_term

    def diff_h(self, state, action, diff_var):
        """
            Derivative of h is either energy if action corresponds to diff_var or zero
        """

        f_num = len(self.f_array.keys())
        if f_num * self.actions.index(action) * 2 <= diff_var < f_num * (self.actions.index(action) + 1) * 2:
            # diff for theta * energy^(1/3)
            if diff_var % 2 == 0:
                features = self.features(state, action)
                if self.h_type == 'linear':
                    return features[int(diff_var / 2)]
                else:
                    return np.sign(features[int(diff_var / 2)]) * (np.abs(features[int(diff_var / 2)])) ** (1 / 3)
            else:
                # diff for theta
                return self.free_term
        else:
            return 0

    def alpha_func(self):
        """
        Returns dynamic learning rate
        """

        x = 5000

        t = (x) / (100 * x + self.step)
        return t

    def load_params(self):
        """
        Loads all parameters from DB
        """

        row = self.sql_execute(self.sql_select_theta, result='one')

        if row is not None:
            self.theta = json.loads(row[0])
        else:
            print(f'Reading theta from DB is failed')

        # load energy params
        if Reinforce.f_energy in self.f_array.keys():
            rows = self.sql_execute(self.sql_select_energy_history + str(self.step_energy_limit),
                                    result='all')
            if rows is not None:
                energies = [item[0] for item in rows]
                if len(energies) > 0:
                    self.f_array[Reinforce.f_energy] = energies
                    if Reinforce.f_old_energy in self.f_array.keys():
                        self.f_array[Reinforce.f_old_energy] = energies
            else:
                print(f'Reading energies from DB is failed')

        self.calc_f_scaling()

        if self.reg_params['scale_reward']:
            if 'step_reward_limit' in self.reg_params:
                limit = self.reg_params['step_reward_limit']
            else:
                limit = self.step_energy_limit
            rows = self.sql_execute(self.sql_select_reward_history + str(limit),
                                    result='all')
            if rows is not None:
                rewards = [item[0] for item in rows]
                if len(rewards) > 1:
                    self.mean_reward = statistics.mean(rewards)
                    self.std_reward = np.std(rewards)
                else:
                    self.mean_reward = None
                    self.std_reward = None
            else:
                self.mean_reward = None
                self.std_reward = None
                print(f'Reading rewards from DB is failed')

    def print_probabilities(self, state, excluded_actions=[]):
        """
        Prints probabilities for actions given the state
        """

        draws = {k: 0 for k in self.actions}

        for k, v in draws.items():
            draws[k] = round(self.action_prob(k, state=state, excluded_actions=excluded_actions), 3)

        print("For energy %f the probabilities are" % state.energy)

    def draw_actions_probabilities(self, min_energy, max_energy, step=0.001,
                                   filename="reinforce_actions_probabilities.png",
                                   title='Actions preferences',
                                   excluded_actions=[],
                                   max_y=100,
                                   small=False):
        """
        Draws the probabilities distribution for actions depending on the energy of the state
        """

        plt.figure(figsize=(800 / 150, 600 / 150), dpi=150)

        if small:
            fontsize = 24
            tickssize = 18
            ratio = 1
            plt.locator_params(nbins=4)
            linewidth = 3
        else:
            fontsize = 13
            tickssize = 13
            ratio = 6 / 8
            plt.locator_params(nbins=6)
            linewidth = 2.5

        x = np.arange(min_energy, max_energy, step)

        plt.rcParams['font.family'] = "sans-serif"
        plt.rcParams['font.sans-serif'] = "Arial"

        for action in self.actions:
            if action not in excluded_actions:
                y = []
                for i in x:
                    state = State(i)
                    prob = self.action_prob(action=action, state=state, excluded_actions=excluded_actions)
                    prob = prob * 100 * (1 - self.reg_params['epsilon']) \
                           + 1 / len(self.actions) * 100 * self.reg_params['epsilon']
                    y.append(prob)

                line, = plt.plot(x, y, linewidth=linewidth)
                line.set_label(action)

        axes = plt.gca()
        axes.set_ylim([0, max_y])
        # if not small:
        plt.legend(prop={'size': 14}, loc='upper right', frameon=False)

        plt.title(title, fontsize=fontsize)
        plt.xlabel('Energy (eV/atom)', fontsize=tickssize)
        plt.ylabel('Probability (%)', fontsize=tickssize)
        plt.xticks(fontsize=tickssize)
        plt.yticks(fontsize=tickssize)

        x_left, x_right = plt.xlim()
        y_low, y_high = plt.ylim()
        plt.gca().set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)

        plt.savefig(filename,
                    bbox_inches='tight',
                    pad_inches=0.1)
        plt.clf()

    def entropy(self, state):
        """
        Returns entropy of the  probabilities distribution on actions
        This will be replace by the built-in function
        """

        res = 0

        for a in self.actions:
            prob = self.action_prob(a, state)
            if prob == 0:
                prob = 0.00001
            try:
                res += prob * math.log(prob)
            except ValueError as e:
                print(e)

        return -res

    def animated_actions_probabilities(self, min_energy, max_energy, filename, interval=10,
                                       sl=[], old_energy=0, max_y=100):
        """
        Shows how action probabilities distribution changes in time
        """

        # if there is no db we work with then we don't need to load theta
        rows = self.sql_execute(self.sql_select_theta_history, result='all')
        thetas = []
        for row in rows:
            thetas.append(json.loads(row[0]))

        if old_energy == 0:
            old_energy = min_energy

        # General plot parameters
        mpl.rcParams['font.family'] = 'Avenir'
        mpl.rcParams['font.size'] = 16
        mpl.rcParams['axes.linewidth'] = 2
        mpl.rcParams['axes.spines.top'] = False
        mpl.rcParams['axes.spines.right'] = False
        mpl.rcParams['xtick.major.size'] = 10
        mpl.rcParams['xtick.major.width'] = 2
        mpl.rcParams['ytick.major.size'] = 10
        mpl.rcParams['ytick.major.width'] = 2
        mpl.rcParams["legend.loc"] = 'center right'

        actions_num = len(self.actions)

        fig, ax = plt.subplots()

        if len(sl) > 0:
            lines = [0] * (actions_num * len(sl))
        else:
            lines = [0] * (actions_num + 1)

        for i in range(len(lines)):
            lines[i], = plt.plot([min_energy, max_energy], [0, max_y])

        T = np.linspace(1, len(thetas), int(len(thetas) / interval))

        temp = ax.text(int((min_energy + max_energy) / 2), 100, '', ha='right', va='top', fontsize=20)

        plt.title("Actions probability distribution")
        plt.xlabel("Energy (eV/atom)")
        plt.ylabel("Probability (%)")

        # colors = plt.get_cmap('gist_ncar', len(lines) + 3)
        colors = plt.get_cmap('tab10', len(lines) + 3)

        remember_theta = self.theta

        if len(sl) > 0:
            # Animation function
            def animate(i):
                x = np.linspace(min_energy, max_energy, 100)
                for local_sl in sl:
                    for a in range(actions_num):
                        y = np.ones(len(x))
                        for j in range(len(x)):
                            self.theta = thetas[int(T[i] - 1)]
                            y[j] = 100 * self.action_prob(self.actions[a],
                                                          State(x[j], np.zeros(local_sl), old_energy=old_energy))
                        lines[len(sl) * a + sl.index(local_sl)].set_data(x, y)
                        lines[len(sl) * a + sl.index(local_sl)].set_label(self.actions[a] + ' sl=' + str(local_sl))
                        lines[len(sl) * a + sl.index(local_sl)].set_color(
                            colors(len(sl) * a + sl.index(local_sl)))

                plt.legend()

                temp.set_text(str(int(T[i])) + ' steps')
                # temp.set_color(colors(i))
        else:
            # Animation function
            def animate(i):
                x = np.linspace(min_energy, max_energy, 100)
                for a in range(actions_num):
                    y = np.ones(len(x))
                    for j in range(len(x)):
                        self.theta = thetas[int(T[i] - 1)]
                        y[j] = 100 * self.action_prob(self.actions[a],
                                                      State(x[j], old_energy=old_energy))
                    lines[a].set_data(x, y)
                    lines[a].set_label(self.actions[a])
                    lines[a].set_color(colors(a))
                # entropy
                y = np.ones(len(x))
                for j in range(len(x)):
                    self.theta = thetas[int(T[i] - 1)]
                    y[j] = 10 * self.entropy(State(x[j], old_energy=old_energy))
                lines[actions_num].set_data(x, y)
                lines[actions_num].set_label('Entropy')
                lines[actions_num].set_color(colors(actions_num))

                plt.legend()

                temp.set_text(str(int(T[i])) + ' steps')
                # temp.set_color(colors(i))

        # Create animation
        ani = FuncAnimation(fig=fig, func=animate, frames=range(len(T)), interval=50, repeat=True)

        # Ensure the entire plot is visible
        fig.tight_layout()

        # Save and show animation
        ani.save(filename, writer=PillowWriter(fps=2))

        self.theta = remember_theta

        # plt.show()

    def fill_params_db(self):
        """
        Creates and fills the theta_params table in DB with predefined theta
        """

        # create DB for MySQL
        if self.reg_params["storage_type"] == "mysql":
            sql_create_db = f'CREATE DATABASE IF NOT EXISTS {self.params_db["database"]}'
            self.sql_execute(sql_create_db, create_db=True)

        # create reinforce table
        separator = ' double, '
        features = separator.join(self.features_to_tablenames()[0]) + ' double'
        sql_create_reinforce_table = "CREATE TABLE IF NOT EXISTS " + self.reinforce_table + """ (
                                            id INT AUTO_INCREMENT PRIMARY KEY,
                                            step int,
                                            old_energy double,
                                            new_energy double,
                                            action varchar(10),
                                            step_size double,
                                            old_struct varchar(200),
                                            new_struct varchar(200),""" \
                                     + features + "); "
        if self.debug:
            print(sql_create_reinforce_table)

        self.sql_execute(sql_create_reinforce_table)

        reward_str = ''
        if self.reg_params['scale_reward']:
            reward_str = ', unscaled_reward double, scaled_reward double'
        # create theta table
        sql_create_theta_table = f"""CREATE TABLE IF NOT EXISTS {self.theta_table} (
                                            id INT AUTO_INCREMENT PRIMARY KEY,
                                            theta text,
                                            reinforce_id tinyint {reward_str},
                                            opt INT
                                        ); """

        self.sql_execute(sql_create_theta_table)

        # if the theta table is empty fill it with default value
        row = self.sql_execute(f'SELECT COUNT(*) FROM {self.theta_table}', result='one')
        if row is None:
            print(f'Reading the number of thetas in DB is failed')
        else:
            if row[0] < 1:
                self.sql_execute(self.sql_update_theta())


    def sql_execute(self, sql, max_attempts=5, result='none', create_db=False):

        if self.reg_params["storage_type"] == "mysql":

            return self.mysql_execute(sql, max_attempts=max_attempts, result=result, create_db=create_db)
        
        else:

            return self.sqlite_execute(sql, max_attempts=max_attempts, result=result, create_db=create_db)
        

    def mysql_execute(self, sql, max_attempts=5, result='none', create_db=False):
        """
        Opens connection to db, executes the query in <max_attempts> attemps, closes connection
        result parameter can be 'one', 'all', if we do the select statement, or 'none' if the query is update
        """

        failed = True
        attempts = 0
        from mysql.connector import connect, Error
        while attempts < max_attempts and failed:
            try:
                # conn_db = sqlite3.connect(self.params_db)
                conn_db = connect(
                    host=self.params_db['host'],
                    user=self.params_db['user'],
                    password=self.params_db['password'],
                    database='' if create_db else self.params_db['database'],
                    auth_plugin='mysql_native_password'
                )
                try:
                    cur = conn_db.cursor(buffered=(result != 'none'))
                    cur.execute(sql)
                    conn_db.commit()
                    failed = False
                    if self.debug:
                        print(f'DB query {sql}: success')

                    if result == 'one':
                        return cur.fetchone()
                    elif result == 'all':
                        return cur.fetchall()
                    cur.close()
                except Error as e:
                    attempts += 1
                    print(f'DB query {sql} failed in attempt {attempts} with error: {e}')
                    time.sleep(10)

                conn_db.close()

            except Error as e:
                attempts += 1
                print(f'DB query {sql} failed in attempt {attempts} with error: {e}')
                time.sleep(10)


    def sqlite_execute(self, sql, max_attempts=5, result='none', create_db=False):
        """
        Opens connection to db, executes the query in <max_attempts> attemps, closes connection
        result parameter can be 'one', 'all', if we do the select statement, or 'none' if the query is update
        """

        failed = True
        attempts = 0

        while attempts < max_attempts and failed:
            try:
                # create connection if 1) database is a file or 2) database is in memory but this is the first call
                if self.params_db['database'] != ":memory:" or self.conn_db == None:
                    self.conn_db = sqlite3.connect(self.params_db['database'])

                try:
                    cur = self.conn_db.cursor()
                    # cur = conn_db.cursor(buffered=(result != 'none'))
                    cur.execute(sql)
                    self.conn_db.commit()
                    failed = False
                    if self.debug:
                        print(f'DB query {sql}: success')

                    if result == 'one':
                        return cur.fetchone()
                    elif result == 'all':
                        return cur.fetchall()
                    cur.close()
                except sqlite3.Error as e:
                    attempts += 1
                    print(f'DB query {sql} failed in attempt {attempts} with error: {e}')
                    time.sleep(10)

                if self.params_db['database'] != ":memory:":
                    self.conn_db.close()

            except sqlite3.Error as e:
                attempts += 1
                print(f'DB query {sql} failed in attempt {attempts} with error: {e}')
                time.sleep(10)
