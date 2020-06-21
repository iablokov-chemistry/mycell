import os
import subprocess
import time
import json
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import sys

# definitions
ps_to_s = 10**-12
A3_to_L = 10**-27
N_a = 6.0 * 10**23

conc_setup = {}
conc_setup["rO2* + LH => L* + rO2H"]       = {"rO2*": 10, "LH": 100}
conc_setup["LO2* + LH => LO2H + L*"]       = {"LO2*": 10, "LH": 100}
conc_setup["LO2* + NOH => NO* + LO2H"]     = {"LO2*": 20, "NOH": 20}
conc_setup["NOH + HO2* => NO* + H2O2"]     = {"NOH": 50, "HO2*": 50}
conc_setup["L* + NO* => LNO"]              = {"NO*": 20, "L*": 20}
conc_setup["HO2* + LH => L* + H2O2"]       = {"HO2*": 10, "LH": 100}
conc_setup["LO2* + NO* => NOH + O2 + L_H"] = {"LO2*": 10, "NO*": 10}
conc_setup["LO2* + LO2* => LO4L"]          = {"LO2*": 20}
conc_setup["HO2* + HO2* => H2O2 + O2"]     = {"HO2*": 40}
conc_setup["NO* + HO2* => NOH + O2"]       = {"HO2*": 30, "NO*": 30}
conc_setup["L* + O2 => LO2*"]              = {"O2": 10, "L*": 10}

conc_setup["LO2* => HO2* + L_H"]           = {"LO2*": 100}


L_medium  = 200
R_micelle = 23
prob_limit = 50
steps = 100000000

n_term_steps = 0

def get_k_1(_df, _r1, _p1, _phase):

    if _phase == "medium":
        # V box
        V = L_medium ** 3               # A^3
    else:
        # V micelle
        V = 4 / 3 * 3.14 * R_micelle ** 3

    #df = _df[_df[_p1] > 0]
    df = _df
    tm = (df.index * ps_to_s).to_numpy().astype(np.float64)

    c_r1 = ((df[_r1] / N_a) / (V * A3_to_L)).to_numpy().astype(np.float64)
    c_p1 = ((df[_p1] / N_a) / (V * A3_to_L)).to_numpy().astype(np.float64)

    x = np.array(tm)
    x = x[:, np.newaxis]
    a, r, _, _ = np.linalg.lstsq(x, c_p1)

    k = a / c_r1
    k_avg = np.round(np.mean(k), 4)
    k_res = np.round((r[0] ** 0.5) / len(tm), 4)

    '''
    p_p1 = np.poly1d(np.polyfit(tm, c_p1, 2))
    p_p1_d = p_p1.deriv()

    k = p_p1_d(tm) / c_r1
    k_avg = np.round(np.mean(k), 4)
    k_std = np.round(np.std(k), 4)
    '''
    return {"k_avg": k_avg, "k_res": k_res}

def get_k_2(_df, _r1, _r2, _p1, _phase):

    if _phase == "medium":
        # V box
        V = L_medium ** 3               # A^3
    else:
        # V micelle
        V = 4 / 3 * 3.14 * R_micelle ** 3

    #df = _df[_df[_p1] > 0]
    df = _df
    tm = (df.index * ps_to_s).to_numpy().astype(np.float64)

    c_r1 = ((df[_r1] / N_a) / (V * A3_to_L)).to_numpy().astype(np.float64)
    c_r2 = ((df[_r2] / N_a) / (V * A3_to_L)).to_numpy().astype(np.float64)
    c_p1 = ((df[_p1] / N_a) / (V * A3_to_L)).to_numpy().astype(np.float64)

    x = np.array(tm)
    x = x[:, np.newaxis]
    a, r, _, _ = np.linalg.lstsq(x, c_p1)

    k = a / (c_r1 * c_r2)
    k_avg = np.round(np.mean(k), 4)
    k_res = np.round((r[0] ** 0.5) / len(tm), 4)

    '''
    p_p1 = np.poly1d(np.polyfit(tm, c_p1, 2))
    p_p1_d = p_p1.deriv()

    k = p_p1_d(tm) / (c_r1 * c_r2)
    k_avg = np.round(np.mean(k), 4)
    k_std = np.round(np.std(k), 4)
    '''
    return {"k_avg": k_avg, "k_res": k_res}


main_config = json.load(open("./cfg/config.main.json", "r"))
rxn_config = json.load(open("./cfg/config.rxn.json", "r"))

k_stat_data = []

for rxn_item in rxn_config["reactions"]:

    if len(sys.argv) == 2:
        if rxn_item["id"] != sys.argv[1]:
            continue

    prob = rxn_item["probability"]

    print(rxn_item["name"])
    k_stat_item = {}
    k_stat_item["rxn"] = rxn_item["name"]
    k_stat_item["k_0"] = rxn_item["k"]

    for t in range(0,3):

        while True:
            cfg = dict(main_config)
            cfg["run"]["steps"] = steps
            #cfg["run"]["steps"] = min(int((1.00E+5 / rxn_item["k"])) * 10, 5000) + 1000

            cfg["run"]["probability-mode"]["enabled"] = True
            cfg["run"]["probability-mode"]["product"] = rxn_item["products"][0]
            cfg["run"]["probability-mode"]["limit"] = prob_limit

            for particle_type in cfg["particle-types"]:
                if rxn_item["phase"] == "medium":
                    particle_type["micelle-bound"] = False
                else:
                    particle_type["micelle-bound"] = True

                if particle_type["type"] in rxn_item["reactants"]:
                    particle_type["fixed-concentration"] = True      # fixed for reactants
                    particle_type["final-product"] = False
                else:
                    particle_type["fixed-concentration"] = False
                    particle_type["final-product"] = True            # to skip brownian dynamics

            if rxn_item["phase"] == "medium":
                cfg["micelle"]["geometry"]["r"] = 0.01
            else:
                cfg["micelle"]["geometry"]["r"] = R_micelle

            cfg["reactions"] = []
            cfg["particles"] = []

            cfg["reactions"].append(dict(rxn_item))
            cfg["reactions"][0]["probability"] = prob

            for particle_type in rxn_item["reactants"]:
                if rxn_item["phase"] == "medium":
                    amount = 20
                    if rxn_item["name"] in conc_setup:
                        amount = conc_setup[rxn_item["name"]][particle_type]
                    cfg["particles"].append({"type": particle_type, "position": "random", "amount": amount})
                else:
                    amount = 20
                    if rxn_item["name"] in conc_setup:
                        amount = conc_setup[rxn_item["name"]][particle_type]
                    cfg["particles"].append({"type": particle_type, "position": "random_micelle", "amount": amount})

            for particle_type in rxn_item["products"]:
                if rxn_item["phase"] == "medium":
                    cfg["particles"].append({"type": particle_type, "position": "random", "amount": 0})
                else:
                    cfg["particles"].append({"type": particle_type, "position": "random_micelle", "amount": 0})

            path_cfg = f"./config/config_{rxn_item['id']}.json"
            json.dump(cfg, open(path_cfg, "w"), indent=4)

            # prepare dir
            cwd = os.getcwd()
            
            os.makedirs("./log/", exist_ok=True)
            threadId = "th_" + rxn_item["id"]
            shutil.rmtree("./log/" + threadId, ignore_errors=True)
            os.makedirs("./log/" + threadId, exist_ok=True)
            os.makedirs("./log/" + threadId + "/coords", exist_ok=True)

            # run micelle
            fLog = open("./log/" + threadId + "/log.txt", "w")
            thread = subprocess.Popen(["./mycell", threadId, path_cfg], cwd=cwd, stdout=fLog)
            thread.wait()
            fLog.close()

            fn_conc = "./log/" + threadId + "/conc.csv"

            # if micelle failed
            if not os.path.exists(fn_conc): continue

            # parse k
            shutil.copy(fn_conc, f"./output/conc/{rxn_item['id']}_{t}.csv")

            df = pd.read_csv(fn_conc, index_col="Time")
            if len(rxn_item["reactants"]) == 1:
                k_data = get_k_1(df, cfg["particles"][0]["type"], cfg["particles"][1]["type"], rxn_item["phase"])
            elif len(rxn_item["reactants"]) == 2:
                k_data = get_k_2(df, cfg["particles"][0]["type"], cfg["particles"][1]["type"], cfg["particles"][2]["type"], rxn_item["phase"])
            else:
                print("Wrong number of reactants")
                exit(0)

            print(cfg["run"]["steps"], prob,  k_data)

            k_stat_item[f"k_avg_{t}"] = round(k_data["k_avg"], 7)
            k_stat_item[f"k_res_{t}"] = round(k_data["k_res"], 7)
            k_stat_item[f"prob_{t}"] = round(prob, 10)

            if k_data["k_avg"] > 0.0:
                prob *= (rxn_item["k"] / k_data["k_avg"])
            else:
                prob *= 10

            break
    print(k_stat_item)
    k_stat_data.append(k_stat_item)

    df_k_stat = pd.DataFrame(k_stat_data)
    df_k_stat.to_csv("./output/prob_stat." + rxn_item["id"] + ".csv")