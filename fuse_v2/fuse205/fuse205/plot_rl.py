import sqlite3
import sys
import math
from rlcsp import Reinforce, State
import numpy as np
import matplotlib.pyplot as plt
import datetime

def standardize_feature_ben(mean_energy, std_energy, val, default_val=0):
    if mean_energy is not None and std_energy != 0 and std_energy is not None:
        return (val - mean_energy) / std_energy
    else:
        return default_val

def features_ben(state, action, actions, mean_energy, std_energy):
    features = ['energy']
    res = np.zeros(len(features) * len(actions))
    f_counter = 0
    res[1 * actions.index(action) + f_counter] = standardize_feature_ben(mean_energy, std_energy,state.energy)
    f_counter += 1
    return res

def h_ben(state, h_type, action, actions, theta, free_term, mean_energy, std_energy):
    features = features_ben(state, action, actions, mean_energy, std_energy)
    #print(f"features: {features}")
    start_coef = []
    for i in range(len(features)):
        if features[i] != 0:
            start_coef.append(theta[2 * i + 1])
            if h_type != 'linear':
                features[i] = np.sign(features[i]) * (np.abs(features[i])) ** (1 / 3)
                
    return np.multiply(theta[::2], features).sum() + sum(start_coef) * free_term

def plot_rl_graph(actions,ini_pop,params_db,alpha,reinforce_table,theta_table,reinforce_id,reg_params,reinforce_verbosity,excluded_actions,h_type,
                          reward_type=Reinforce.change_in_features,
                          features_set=['energy'],
                          episode_length=1,
                          max_energy=0,):

    #get max and min energy from .db file and means and std
    mean_energy_list=[]
    std_energy_list=[]
    
    #print(actions)
    #sys.exit()
    
    # if sqlite is being used
    if reg_params['storage_type'] == 'sqlite':
        conn = sqlite3.connect('../'+str(params_db['database']))
        cursor = conn.cursor()
        try:
            cursor.execute(f'''
                  SELECT max_energy
                  FROM {reinforce_table}
                   ORDER BY id DESC
            LIMIT 1
            ''')
            result = cursor.fetchone()
            if result:
                max_energy=float(result[0])
            else:
                print('max energy not found')
            cursor.execute(f'''
                  SELECT min_energy
                  FROM {reinforce_table}
                   ORDER BY id DESC
            LIMIT 1
            ''')
            result = cursor.fetchone()
            if result:
                min_energy =float(result[0])
            else:
                print('min energy not found')
            cursor.execute(f'''
                      SELECT mean_energy
                      FROM {reinforce_table}
                       ORDER BY id ASC
                    ''')
            rows = cursor.fetchall()
            for row in rows:
                mean_energy = row[0]
                mean_energy_list.append(mean_energy)
            cursor.execute(f'''
                              SELECT std_energy
                              FROM {reinforce_table}
                               ORDER BY id ASC
                            ''')
            rows = cursor.fetchall()
            for row in rows:
                std_energy = row[0]
                std_energy_list.append(std_energy)
        finally:
            conn.close()
        theta_list=[]
        conn = sqlite3.connect('../'+str(params_db['database']))
        cursor = conn.cursor()
        try:
            cursor.execute(f'''
              SELECT theta
              FROM {theta_table}
                WHERE opt = 1
                ORDER BY id ASC
            ''')
            rows = cursor.fetchall()
            for row in rows:
                theta=row[0]
                theta = theta.strip('[]')
                theta = [float(x) for x in theta.split(',') if x]
                theta_list.append(theta)
            cursor.execute(f"SELECT COUNT(*) FROM {theta_table};")
        except sqlite3.Error as e:
            print(f"An error occurred: {e}")
        finally:
            conn.close()

    if reg_params['storage_type'] == 'mysql':
        from mysql.connector import connect, Error
        conn = connect(host=params_db['host'],
                    user=params_db['user'],
                    password=params_db['password'],
                    database=params_db['database'],
                    auth_plugin='mysql_native_password'
                )
        cursor = conn.cursor()
        try:
            cursor.execute(f'''
                          SELECT max_energy
                          FROM {reinforce_table}
                           ORDER BY id DESC
                    LIMIT 1
                    ''')
            result = cursor.fetchone()
            if result:
                max_energy = float(result[0])
            else:
                print('max energy not found')
            cursor.execute(f'''
                          SELECT min_energy
                          FROM {reinforce_table}
                          ORDER BY id DESC
                     LIMIT 1
                    ''')
            result = cursor.fetchone()
            if result:
                min_energy = float(result[0])
            else:
                print('min energy not found')
            cursor.execute(f'''
                              SELECT mean_energy
                              FROM {reinforce_table}
                               ORDER BY id ASC
                            ''')
            rows = cursor.fetchall()
            for row in rows:
                mean_energy = row[0]
                mean_energy_list.append(mean_energy)
            cursor.execute(f'''
                                      SELECT std_energy
                                      FROM {reinforce_table}
                                       ORDER BY id ASC
                                    ''')
            rows = cursor.fetchall()
            for row in rows:
                std_energy = row[0]
                std_energy_list.append(std_energy)
        finally:
            conn.close()
        theta_list = []
        conn = connect(host=params_db['host'],
                    user=params_db['user'],
                    password=params_db['password'],
                    database=params_db['database'],
                    auth_plugin='mysql_native_password'
                )
        cursor = conn.cursor()
        try:
            cursor.execute(f'''
                      SELECT theta
                      FROM {theta_table}
                        WHERE opt = 1
                        ORDER BY id ASC
                    ''')
            rows = cursor.fetchall()
            for row in rows:
                theta = row[0]
                theta = theta.strip('[]')
                theta = [float(x) for x in theta.split(',') if x]
                theta_list.append(theta)
            cursor.execute(f"SELECT COUNT(*) FROM {theta_table};")
        except sqlite3.Error as e:
            print(f"An error occurred: {e}")
        finally:
            conn.close()

    steps = np.linspace(min_energy,max_energy,50)
    step=steps[1]-steps[0]

    all_probabilities={}
    #used=[]
    for cur_action in actions:
        #print(f"\n\naction: {cur_action}")
        
        #if cur_action in used:
        #	  continue
        	  
        #used.append(cur_action)
        
        probabilities_2 = {}
        count = 0
        for i in theta_list:
            count+=1
            h_dict={str(key): [] for key in actions}
            #print(f"h_dict: {h_dict}")
            #sys.exit()
            probs=[]
            y = np.arange(min_energy, max_energy, step)
            #print(f"y: {y}")
            used=[]
            for k in actions:                
                if k in used:
                	 continue
                used.append(k)
                #print(f"k: {k}")
                for j in y:
                    #print(f"j: {j}")
                    state = State(j)
                    #print(f"state: {state}")
                    #print(f"h_type: {h_type}")
                    #print(f"actions: {actions}")
                    #print(f"theta: {theta}")
                    #print(f"mean energy: {mean_energy_list[(count-2)+ini_pop]}")
                    #print(f"std_energy: {std_energy_list[(count-2)+ini_pop]}")
                    
                    h = h_ben(state=state, h_type=h_type, action=k, actions=actions, theta=i, free_term=1, mean_energy=mean_energy_list[(count-2)+ini_pop], std_energy=std_energy_list[(count-2)+ini_pop])
                    h_dict[str(k)].append(h)

            ## need a prob for each action, given the hs of the other actions in that state
            #print(h_dict)
            #print(h_dict['1'])
            
            #for k in list(h_dict.keys()):
            #	print(len(h_dict[k]))
            
            #sys.exit()
            #used=[]
            for k in range(len(h_dict['1'])):
                beta = 1
                
                for a in actions:
                    if a in excluded_actions:
                        return 0
                       
                h_list = []
                h = 0
                #print(actions)
                #print(h_dict.keys())
                
                for a in actions:
                    if a not in excluded_actions:
                        h_list.append(beta * h_dict[str(a)][k])
                        if a == cur_action: ####### this bit selects the action for the plot
                            h = h_list[-1]
                            
                h_list = np.array(h_list)
                max_h = h_list.max()
                h -= max_h
                h_list -= max_h
                prob=math.exp(h) / np.exp(h_list).sum()
                probs.append(prob * 100)

            probabilities_2[count] = probs
        all_probabilities[str(cur_action)] = probabilities_2
       
    #print(all_probabilities.keys())  
    min_prob=100
    max_prob=0
    for i in all_probabilities.keys():
        for k in all_probabilities[i].keys():
            if max(list(all_probabilities[i][k])) > max_prob:
                max_prob=max(list(all_probabilities[i][k]))
            
            if min(list(all_probabilities[i][k])) < min_prob:
                min_prob=min(list(all_probabilities[i][k]))

    for i in all_probabilities.keys():
        x_vals = []
        y_vals = []
        color_vals = []
        for k in all_probabilities[i].keys():
            x_vals.extend([k] * len(y))  # Repeating 'i' for each y-value
            y_vals.extend(y)  # Adding all y-values for this 'i'
            color_vals.extend(all_probabilities[i][k])  # Adding the corresponding color values

            # Now scatter plot all at once
        plt.style.use('fast')
        plt.scatter(x_vals, y_vals, c=color_vals, cmap='jet', s=20, alpha=0.5)
        plt.colorbar(label='Probability')
        plt.clim(min_prob,max_prob)
        tic=int(count/10)
        plt.xticks(np.arange(1, count, step=tic))
        plt.ylabel('Energy eV/atom')
        plt.xlabel('Steps')
        plt.savefig('probability'+str(i)+'.png',dpi=600)
        plt.clf()
        plt.close()
